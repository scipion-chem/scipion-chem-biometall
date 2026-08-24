#!/usr/bin/python3
"""
homology_search.py

Utility module for homology-based metal identification (Path 1).

Queries the RCSB PDB for sequence-similar structures that contain a metal
ion of interest.  Returns the metal type found in the best homologue so
that MetalKB can place it in the query structure.

No structural alignment or coordinate transfer is performed here: MetalKB
computes the optimal ion position directly using MESPEUS statistical
potentials for the identified metal, which is more accurate than transferring
coordinates from a different protein context.

Public API

search_homologs(sequence, ...)  ? list[dict]  ranked homologue hits
"""

import json
import urllib.request
import urllib.error

# Metal set
METALS_OF_INTEREST = {'ZN', 'CA', 'MG', 'FE', 'CU', 'MN', 'NI'}


RCSB_SEARCH_URL  = 'https://search.rcsb.org/rcsbsearch/v2/query'
RCSB_GRAPHQL_URL = 'https://data.rcsb.org/graphql'

# Sequence-based homologue search

def search_homologs(sequence, identityCutoff=0.30, evalueCutoff=1e-5,
                    maxResults=10, metals=None, timeout=30):
    """
    Query the RCSB PDB REST API for structures homologous to *sequence* that
    contain at least one metal ion of interest.

    Uses the RCSB sequence-similarity service (MMseqs2 internally), then
    filters hits via GraphQL to confirm the presence of the required metal.

    Parameters

    sequence        : str   Amino-acid sequence (no FASTA header).
    identity_cutoff : float Minimum fractional identity (0.30 = 30 %).
    evalue_cutoff   : float Maximum alignment e-value (1e-5 recommended).
    max_results     : int   Maximum number of hits to return.
    metals          : set   Metal symbols to accept.  None ? METALS_OF_INTEREST.
    timeout         : int   HTTP request timeout in seconds.

    Returns

    list[dict]
        Sorted by identity (descending).  Each entry:
          'pdb_id'   : str    ? 4-letter PDB accession code
          'chain_id' : str    ? matching polymer chain identifier
          'identity' : float  ? fractional sequence identity (0?1)
          'evalue'   : float  ? alignment e-value (may be None)
          'metals'   : list   ? metal symbols found in this structure
    """
    if metals is None:
        metals = METALS_OF_INTEREST

    # Over-query (max_results × 3) to compensate for post-filter metal losses.
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "sequence",
                    "parameters": {
                        "evalue_cutoff":   evalueCutoff,
                        "identity_cutoff": identityCutoff,
                        "target":          "pdb_protein_sequence",
                        "value":           sequence,
                    }
                }
            ]
        },
        "return_type": "polymer_entity",
        "request_options": {
            "paginate": {"start": 0, "rows": maxResults * 3},
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": "score", "direction": "desc"}],
            "scoring_strategy": "sequence",
        }
    }

    data = json.dumps(query).encode('utf-8')
    req  = urllib.request.Request(
        RCSB_SEARCH_URL,
        data=data,
        headers={'Content-Type': 'application/json'},
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            result = json.loads(resp.read().decode('utf-8'))
    except urllib.error.HTTPError as e:
        raise RuntimeError(f'RCSB search HTTP error {e.code}: {e.reason}')
    except urllib.error.URLError as e:
        raise RuntimeError(f'RCSB search connection error: {e.reason}')

    hits = result.get('result_set', [])
    if not hits:
        return []

    candidates = []
    for hit in hits:
        ident_str = hit.get('identifier', '')
        parts     = ident_str.split('_')
        if len(parts) < 2:
            continue
        pdb_id   = parts[0].upper()
        identity = hit.get('score', 0.0)

        # Extract e-value and chain from nested services/nodes structure
        evalue   = None
        chain_id = 'A'
        for service in hit.get('services', []):
            for node in service.get('nodes', []):
                mp = node.get('match_context', [{}])[0]
                if 'evalue' in mp:
                    evalue = mp['evalue']
                if 'query_entity_mapping' in mp:
                    chains = mp['query_entity_mapping'].get('entity_chains', [])
                    if chains:
                        chain_id = chains[0]

        candidates.append({
            'pdb_id':   pdb_id,
            'chain_id': chain_id,
            'identity': identity,
            'evalue':   evalue,
            'metals':   [],
        })

    filtered = _filter_by_metal(candidates, metals, timeout=timeout)
    filtered.sort(key=lambda h: h['identity'], reverse=True)
    return filtered[:maxResults]


# Internal: metal content filter

def _filter_by_metal(candidates, metals, timeout=30):
    """
    For each candidate PDB, query the RCSB GraphQL API to check whether it
    contains any of the requested metal ions.  Structures without a matching
    metal are discarded.

    Sets c['metals'] = [metal_symbol, ...] (types present)
    and  c['metal_counts'] = {metal_symbol: count} (how many of each).
    """
    valid = []
    for c in candidates:
        try:
            metals_in_pdb = _get_metals_in_structure(c['pdb_id'], timeout=timeout)
        except Exception:
            metals_in_pdb = {}
        found = {m: n for m, n in metals_in_pdb.items() if m in metals}
        if found:
            c['metals']       = list(found.keys())
            c['metal_counts'] = found
            valid.append(c)
    return valid


def _get_metals_in_structure(pdb_id, timeout=30):
    """
    Query RCSB GraphQL for the non-polymer components in *pdb_id*.
    Returns a dict {metal_symbol: count} for metals in METALS_OF_INTEREST.
    Example: {'ZN': 4} for alcohol dehydrogenase (1CDO).
    """
    graphql_query = """
    {
      entry(entry_id: "%s") {
        nonpolymer_entities {
          rcsb_nonpolymer_entity {
            pdbx_number_of_molecules
          }
          nonpolymer_comp {
            chem_comp {
              id
              type
            }
          }
        }
      }
    }
    """ % pdb_id.upper()

    data = json.dumps({'query': graphql_query}).encode('utf-8')
    req  = urllib.request.Request(
        RCSB_GRAPHQL_URL,
        data=data,
        headers={'Content-Type': 'application/json'},
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            result = json.loads(resp.read().decode('utf-8'))
    except Exception:
        return {}

    metals_found = {}
    try:
        entities = (result.get('data', {})
                          .get('entry', {})
                          .get('nonpolymer_entities') or [])
        for ent in entities:
            comp    = (ent.get('nonpolymer_comp') or {})
            cc      = (comp.get('chem_comp') or {})
            chem_id = (cc.get('id') or '').upper()
            if chem_id in METALS_OF_INTEREST:
                count = ((ent.get('rcsb_nonpolymer_entity') or {})
                             .get('pdbx_number_of_molecules') or 1)
                metals_found[chem_id] = int(count)
    except Exception:
        pass

    return metals_found