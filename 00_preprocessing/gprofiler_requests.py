import requests
import pandas as pd

def request_gost_link(members):
    r = requests.post(
        url='https://biit.cs.ut.ee/gplink/l',
        json={
        'url': 'https://biit.cs.ut.ee/gprofiler/gost',
        'payload': {
            'organism':'mmusculus',
            'query': members,
            'sources' :["GO:MF","GO:CC","GO:BP","KEGG","REAC","WP","TF","MIRNA","HPA","CORUM","HP"], 
            'no_evidences':False, #skip lookup for evidence codes. Speeds up queries, if there is no interest in evidence codes.
            'no_iea':False, #Ignore electonically annotated GO annotations,
            'ordered': False, #not an ordered query
        },
#         'data_version': "e106_eg53_p16_65fcd97"
        
        },
        headers={
        'User-Agent':'FullPythonRequest'
        }
    )
    return r


def request_gost(members):
    r = requests.post(
        url='https://biit.cs.ut.ee/gprofiler/api/gost/profile/',
        json={
        'organism':'mmusculus',
        'query': members,
        'sources' :["GO:MF","GO:CC","GO:BP","KEGG","REAC","WP","TF","MIRNA","HPA","CORUM","HP"], 
        'no_evidences':False, #skip lookup for evidence codes. Speeds up queries, if there is no interest in evidence codes.
        'no_iea':False, #Ignore electonically annotated GO annotations,
        'ordered': True #not an ordered query
        },
        headers={
        'User-Agent':'FullPythonRequest'
        }
    )
    return r


def make_gost_df(module_members):
    module_name, members = [x for x in module_members.items()][0]
    cols = ['module',
            'link',
             'native',
             'source',
             'name',
             'description',
             'p_value',
             'significant',
             'term_size',
             'all_module_genes',
             'intersections',
             'intersection_size',
             'query_size',       
             'effective_domain_size',
             'source_order',
             'goshv',
             'parents',
             'precision',
             'recall',
             'query',
             'group_id',]

    gost_df_cols = ['description',
                    'effective_domain_size',
                    'goshv',
                    'intersection_size',
                    'intersections',
                    'name',
                    'native',
                    'p_value',
                    'parents',
                    'precision',
                    'query',
                    'query_size',
                    'recall',
                    'significant',
                    'source',
                    'term_size',
                    'source_order',
                    'group_id']

    r = request_gost(members)

    #for module_name, r in gost_requests.items():
    #     if 'miss' in module_name:
    #         continue
    gost_df = pd.DataFrame(r.json()['result'])
    if len(gost_df) < 1:
        gost_df = pd.DataFrame({x:['NA'] for x in gost_df_cols})
    gost_df['module'] = module_name
    if 'intersections' in gost_df.columns.tolist():
        intersections = gost_df['intersections'].tolist()
        intersections_remapped = []
        all_genes = []
        for ix in intersections:
            ix = [x if x is not None else [] for x in ix]
            ix_idx = [i for i, e in enumerate(ix) if len(e) > 0]
            members = [module_members[module_name][x] for x in ix_idx]
            intersections_remapped.append(members)
            all_genes.append(module_members[module_name])
        gost_df['intersections'] = intersections_remapped
        gost_df['all_module_genes'] = all_genes
    gost_df['link'] = 'https://biit.cs.ut.ee/gplink/l/'+ request_gost_link(members).json()['result']
    gost_df = gost_df[cols]
    return gost_df