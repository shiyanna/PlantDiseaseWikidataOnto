#%%
# Source code for "Automatic Information Retrieval and Extraction Methodology for the Ontology of Plant Diseases"

#%%
import sys
import os
from typing import Annotated, Dict, List, Literal, Optional, Tuple, Union
import requests
from SPARQLWrapper import SPARQLWrapper, JSON
import json
import pandas as pd
import numpy as np
from tqdm import tqdm
import time
from functools import lru_cache
from bs4 import BeautifulSoup
from fuzzywuzzy import process
from json_repair import repair_json

#%%
endpoint_url = "https://query.wikidata.org/sparql"

#%%
t_q = str
t_q_guaranteed = str
t_eppo = str
t_eppo_guaranteed = str
t_sparql_response = List[Optional[Dict]]
t_raw_sparql_response = Dict[str,Dict[str,List]]
t_triplet = Tuple[str,Optional[t_q],Optional[t_eppo]]
t_triplet_list = List[t_triplet]
#%%
# tuplify = lambda arr: tuple(i if type(i) == str else tuple(i) for i in arr if i) if arr else ()

tuplify = lambda arr: tuple(i if type(i) in (str,type(None)) else tuple(i) for i in arr if i or i is None) if arr else ()
#%%
#=========
# EPPO


#%%

# EPPO api token. Get it here: https://data.eppo.int
_AUTHTOKEN = 'your-token-here'

#%%
# 0/0

# Get downloadable eppo files here: https://data.eppo.int/
# You need one called Bayer

@lru_cache(maxsize=None)
def read_local_eppo():
    _df1 = pd.read_csv("eppo/gafname.txt")
    _df2 = pd.read_csv("eppo/gainame.txt")
    _df3 = pd.read_csv("eppo/pflname.txt")
    _df = pd.concat((_df1,_df2,_df3))    
    _names = _df["fullname"].dropna().tolist()
    return _df, _names

@lru_cache(maxsize=None)
def _search_local_eppo_names(name: str, limit) -> List[Tuple[str,int]]:
    _df,_names = read_local_eppo()
    column_as_list = _names
    res = tuplify(process.extract(name, column_as_list,limit=limit))
    
    return res

#%%

def search_local_eppo_name(name: str, limit = 20, _dfindings={}, threshold=90) -> Optional[t_eppo]:
    _eppo_df,_names = read_local_eppo()
    
    if name in _dfindings:
        _res = _dfindings[name]
    else:
        _res = _search_local_eppo_names(name, limit)
    trg = _res[0]
    
    if trg[1] <= threshold:
        return None
    else:
        eppo = _eppo_df[_eppo_df.fullname==trg[0]].code.iloc[0]
        return eppo
        

#%%
def json_post_request(data, url, headers = None) -> dict:
    
    if not headers:
        headers = {
            'User-Agent': "WDQS-example Python/%s.%s" % (sys.version_info[0], sys.version_info[1]),  
            'X-Requested-With': 'XMLHttpRequest',
            'Origin': 'https://data.eppo.int',
            'Connection': 'keep-alive',
        }
    
    response = requests.post(url, headers=headers, data=data)

    if response.status_code == 200:
        try:
            return response.json()
        except ValueError:
            print("Response content is not in JSON format:", response.text)
    else:
        print(f"Request failed with status code {response.status_code}: {response.text}")
        
def json_get_request(url, params, headers = None) -> dict:
    
    if not headers:
        headers = {
            'User-Agent': "WDQS-example Python/%s.%s" % (sys.version_info[0], sys.version_info[1]),  
            'X-Requested-With': 'XMLHttpRequest',
            'Origin': 'https://data.eppo.int',
            'Connection': 'keep-alive',
        }
    
    response = requests.get(url, headers=headers, params=params)

    if response.status_code == 200:
        try:
            return response.json()
        except ValueError:
            print("Response content is not in JSON format:", response.text)
    else:
        print(f"Request failed with status code {response.status_code}: {response.text}")

# Search list of taxons in eppo by api
def _search_taxons(taxons:list,authtoken = None) -> List[Optional[t_eppo]]:
    _url = 'https://data.eppo.int/api/rest/1.0/tools/names2codes'
    if not authtoken:
        global _AUTHTOKEN
        authtoken = _AUTHTOKEN

    if type(taxons) == list:
        intext = "|".join(taxons)
    else:
        intext = str(taxons)
        
    data = {
        'authtoken': authtoken,
        'intext': intext
    }
    
    res_raw = json_post_request(data=data, url=_url)
    s = res_raw["response"]
    res = [i.split(";")[1] for i in s.split("|")]
    res = [i if not "*" in i else None for i in res]
    return res

def search_taxons(taxons:List[str], disable_tqdm = False) -> List[Optional[t_eppo]]:
    # eppo API request
    _res1: List[Optional[t_eppo]] = _search_taxons(taxons) 
    # fuzzy search
    _res2 = [search_local_eppo_name(taxon,_dfindings = _dfindings) for taxon in tqdm(taxons,disable=disable_tqdm)]
    # prioritize fuzzy search
    _res = [i2 or i1 for i1,i2 in zip(_res1,_res2)]
    return _res
    

#%%
0/0 # Be careful not to reexecute to not lose cache

# Get taxon's hosts from eppo by eppo code
@lru_cache(maxsize=None)
def get_host(taxon,authtoken = None) -> Dict[Literal['Major host', 'Host'],List[Dict]]:
    if not authtoken:
        global _AUTHTOKEN
        authtoken = _AUTHTOKEN
        
    params = {
        'authtoken': authtoken
    }

    url = 'https://data.eppo.int/api/rest/1.0/taxon/%s/hosts'
    
    res = json_get_request(url=url%taxon, params=params)
    return res

@lru_cache(maxsize=None)
def get_taxon_info(taxon,authtoken = None):
    if not authtoken:
        global _AUTHTOKEN
        authtoken = _AUTHTOKEN
        
    params = {
        'authtoken': authtoken
    }

    url = 'https://data.eppo.int/api/rest/1.0/taxon/%s'
    
    res = json_get_request(url=url%taxon, params=params)
    return res

@lru_cache(maxsize=None)
def get_taxonomy(taxon,authtoken = None):
    if not authtoken:
        global _AUTHTOKEN
        authtoken = _AUTHTOKEN
        
    params = {
        'authtoken': authtoken
    }

    url = 'https://data.eppo.int/api/rest/1.0/taxon/%s/taxonomy'
    
    res = json_get_request(url=url%taxon, params=params)
    return res
#%%


#%%


#%%

# END EPPO
#=========

#%%

#%%
0/0 # Be careful not to reexecute to not lose cache

# Execute sparql query in wikidata
@lru_cache(maxsize=None)
def get_results(query, endpoint_url = None) -> t_raw_sparql_response:
    if not endpoint_url:
        endpoint_url = "https://query.wikidata.org/sparql"
    user_agent = "WDQS-example Python/%s.%s" % (sys.version_info[0], sys.version_info[1])
    sparql = SPARQLWrapper(endpoint_url, agent=user_agent)
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    return sparql.query().convert()

#%%

#@
### Первый пакет данных (taxons 4)

# Берётся полностью с wikidata
# Сначала берётся класс plant disease (Q2662845) и все его подклассы
# Получаем все (наверное) классы болезней. И много мусора, но это на данном этапе не важно, (потому что у мусора не будет инстансов болезней)

#%%
q_subclasses_of_desease = '''
SELECT ?item ?itemLabel WHERE {
  ?item wdt:P279 wd:Q2662845 .
  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
'''

res = get_results(q_subclasses_of_desease)

subdeseases = (res["results"]["bindings"])
subdeseases_tuples = [[i["itemLabel"]["value"],i["item"]["value"].split("/")[-1]] for i in subdeseases]
subdeseases_qs:List[t_q] = [i[1] for i in subdeseases_tuples]
subdeseases_qs.append("Q2662845")
len(subdeseases_qs)
#%%
#%%

#@
# Далее запрашиваем все instance of (P31) болезней ?DES, где ?DES принимает значения классов из списка `subdeseases_qs`

q_all_deseases = '''
SELECT ?item ?itemLabel WHERE {
  ?item wdt:P31 ?DES .
  VALUES ?DES { %s }.
  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
'''

def f_q_all_deseases(qs): 
    if type(qs) == list:
        s_qs = " ".join([f'wd:{q}' for q in qs])
    else:
        s_qs = f'wd:{qs}'
    res = get_results(q_all_deseases%s_qs)
    return res


#%%

all_deseases: t_sparql_response

#%%
res = f_q_all_deseases(subdeseases_qs)
all_deseases = res["results"]["bindings"]
len(all_deseases)
#%%
#@
# Получаем список болезней. Но и много мусора, который тоже потом отфильтруется (потому что у мусора не будет cause-ов)
#%%

os.makedirs("data", exist_ok=True)

#%%

# checkpoint

with open("data/all_deseases.json", "w+") as fd:
    json.dump(all_deseases, fd, indent=2)

#%%

# with open("data/all_deseases.json", "r+") as fd:
#     all_deseases = json.load(fd)


#%%

all_deseases

#%%
all_deseases_flat:List[Tuple[str, t_q]] = [(i["itemLabel"]["value"],i["item"]["value"].split("/")[-1]) for i in all_deseases]
all_deseases_flat

#%%

#@ 
# Далее получаем список таксонов-патогенов, которые связанны с болезнями (переменная ?DES) 
# отношениями (переменная ?PC): has cause (P828), has immediate cause (P1478) has biological vector (P11231)

#%%
q_all_causes = '''
SELECT ?item ?itemLabel ?caused ?causedLabel WHERE {
  ?item wdt:P31 ?DES ;
        ?PC ?caused .
  VALUES ?DES { %s }.
  VALUES ?PC { wdt:P828 wdt:P11231 wdt:P1478 }.
  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
'''

def f_q_all_causes(qs): 
    if type(qs) == list:
        s_qs = " ".join([f'wd:{q}' for q in qs])
    else:
        s_qs = f'wd:{qs}'
    res = get_results(q_all_causes%s_qs)
    return res

causes = f_q_all_causes(subdeseases_qs)
causes:t_sparql_response = causes["results"]["bindings"]
len(causes)

#%%
causes
#%%

with open("data/causes.json", "w+") as fd:
    json.dump(causes, fd, indent=2)
#%%
# with open("data/causes.json", "r+") as fd:
#     causes = json.load(fd)

#%%

causes_flat:List[Tuple[t_q, str]] = []

for item in causes:
    q = item["caused"]["value"].split("/")[-1]
    label = item["causedLabel"]["value"]
    causes_flat.append((q,label))

deseases_causes_flat: List[Tuple[t_q, str, t_q, str]] = []

for item in causes:
    desease_label = item["itemLabel"]["value"]
    desease_q = item["item"]["value"].split("/")[-1]
    q = item["caused"]["value"].split("/")[-1]
    label = item["causedLabel"]["value"]
    deseases_causes_flat.append((desease_q,desease_label,q,label))

#%%
df_deseases = pd.DataFrame(data=deseases_causes_flat, columns=["desease_q","desease_label","cause_q","cause_label"])
df_deseases
#%%
#@
# Получили список таксонов-патогенов
# Ищем для них eppo коды

# Это делается двумя путями

# Первый быстрый путь - обратиться к API eppo: https://data.eppo.int/api/rest/1.0/tools/names2codes
# Это реализовано в функции _search_taxons. 
# Данный способ иногда даёт ложные коды, вместо точных совпадений выдавая дочерние таксоны, или таксоны имеющие искомую строку в качестве синонима.

# Второй более надёжный но долгий путь - реализовать поиск поверх базы eppo самостоятельно. 
# Это реализовано в функции search_local_eppo_name

# Итоговая функция search_taxons использует оба метода, отдавая предпочтение поиску по базе (второму способу)

#%%
causes_dict:Dict[str,t_q] = {label:q for q,label in causes_flat}
names = list(causes_dict)

#%%
taxons = search_taxons(names)

#%%
causes_triples = [(names[i],causes_dict[names[i]], taxons[i]) for i in range(len(causes_dict))]
causes_triples

df = pd.DataFrame(data = causes_triples, columns=["name", "q", "eppo"])
df
#%%
df.to_csv("data/taxons1.csv",index=False)

#%%
# df = pd.read_csv("data/taxons1.csv")
# df.loc[df[df.eppo.isna()].index,"eppo"] = None

#%%

# По eppo кодам достаём хосты

# За это отвечает этот эндпоинт
# https://data.eppo.int/api/rest/1.0/taxon/%s/hosts
# И функция get_host

# По get запросу с eppo кодом на месте %s он отдаёт json с хостами

#%%

# Makes api call, so can be interrupted on rate limit

hosts:List[Optional[Dict[Literal['Major host', 'Host'],List[Dict]]]] = []
for eppo in tqdm(df["eppo"]):
    hits = get_host.cache_info().hits

    if not eppo:
        hosts.append(None)
        continue
    
    host = get_host(eppo)
    hosts.append(host)
    
    if hits == get_host.cache_info().hits:
        time.sleep(2)
    

hosts
#%%

hosts_tuples: List[Optional[Tuple[t_eppo, str]]] = []

for host in hosts:
    host:Optional[Dict[Literal['Major host', 'Host'],List[Dict]]]
    if not host:
        res = []
    else:
        tuples = []
        for host_item in host.values(): # ['Major host', 'Host']
            for host_i in host_item:
                tuples.append([host_i["eppocode"],host_i["full_name"]])
        res = tuples
    hosts_tuples.append(res)
hosts_tuples

#%%
df['hosts'] = hosts_tuples

#%%

#%%



#%%

#%%
_df_save = df.copy()
_df_save.hosts = _df_save.hosts.apply(json.dumps)
_df_save.to_csv("data/taxons2.csv",index=False)
#%%
# df = pd.read_csv("data/taxons2.csv")
# df.loc[df[df.eppo.isna()].index,"eppo"] = None
# df.loc[df[df.hosts.isna()].index,"hosts"] = '[]'
# df.hosts = df.hosts.apply(json.loads)

#%%

# Вначале мы вытаскивали по болезням патогены . Теперь ставим патогенам в соответствие болезни

corr_deseases: List[Tuple[str, t_q_guaranteed]] = []
for q in tqdm(df["q"]):
    _df = df_deseases[df_deseases["cause_q"]==q][["desease_q","desease_label"]]
    corr_desease = [(a,b) for a,b in zip(_df["desease_label"],_df["desease_q"])]
    corr_deseases.append(corr_desease)
corr_deseases

df["corr_deseases"] = corr_deseases
df

#%%
#@
# Для каждого хоста из eppo нужно найти его q в Wikidata
# Это делается sparql запросом, который ищет таксон по совпадению label-ов
#%%
# q_search_by_label = '''
# SELECT distinct ?item ?itemLabel WHERE{  
#   ?item ?label ?LAB .  
#   ?article schema:about ?item .
#   VALUES ?LAB { %s }
#   SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }    
# }
# '''

# def f_q_search_by_labels(labels): 
#     if type(labels) == list:
#         s_labels = " ".join([f'"{label}"@en' for label in labels])
#     else:
#         s_labels = f'"{labels}"@en'
#     # print(q_search_by_label%s_labels)
#     res = get_results(q_search_by_label%s_labels)
#     return res
get_results
#%%

# # p31 instance of  q16521 taxon

q_search_taxon_by_label = '''
SELECT distinct ?item ?itemLabel WHERE{  
  ?item ?label ?LAB .  
  ?item wdt:P31 wd:Q16521 . 
  ?article schema:about ?item .
  VALUES ?LAB { %s }
  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }    
}
'''

def f_q_search_taxons_by_labels(labels): 
    if type(labels) == list:
        s_labels = " ".join([f'"{label}"@en' for label in labels])
    else:
        s_labels = f'"{labels}"@en'
    # print(q_search_by_label%s_labels)
    res = get_results(q_search_taxon_by_label%s_labels)
    return res

f_q_search_by_labels = f_q_search_taxons_by_labels
#%%
0/0
@lru_cache(maxsize=None)
def f_q_search_taxon_by_label(label:str) -> t_sparql_response: 
    s_labels = f'"{label}"@en'
    res = get_results(q_search_taxon_by_label%s_labels)["results"]["bindings"]
    return res

#%%

def p(items:str, size:int, page:int): # paginator
    return items[page*size:(page+1)*size]
#%%

#%%

hosts:pd.Series = df['hosts']
hosts:List[Optional[List[Tuple[t_eppo,str]]]]

hosts_flat:List[Tuple[t_eppo,str]] = []

for host in hosts:
    if host:
        hosts_flat.extend(host)

#%%

hosts_dict = {k:v for k,v in hosts_flat}
hosts_tuples = list(hosts_dict.items())
hosts_dict_names = {v:k for k,v in hosts_dict.items()}
hosts_names = [i[1] for i in hosts_tuples]

#%%
hosts_qs:List[Optional[t_q]] = []


for host_name in tqdm(hosts_names):
    hits = f_q_search_taxon_by_label.cache_info().hits
    
    res = f_q_search_taxon_by_label(host_name)
    q: t_q = res[0]["item"]["value"].split('/')[-1] if res else None
    hosts_qs.append(q)
    
    if hits == f_q_search_taxon_by_label.cache_info().hits:
        time.sleep(1.4)
 
len(hosts_qs)

#%%

wd_items_dict = {hosts_names[i]:hosts_qs[i] for i in range(len(hosts_names))}
wd_items_dict

#%%

hosts_dict_complete = {name:(hosts_dict_names[name], wd_items_dict[name]) for name in hosts_dict_names}
hosts_dict_complete

#%%

# Form new hosts column for dataframe
hosts_col:List[List[Tuple[str,t_q,t_eppo]]] = []
for host in df["hosts"]:
    host:List[Tuple[t_eppo,str]]
    if not host:
        hosts_col.append([])
        continue
    new_host = []
    for host_i in host:
        name = host_i[1]
        eppo = host_i[0]
        q = hosts_dict_complete[name][1]
        triple = (name,q,eppo)
        new_host.append(triple)
    hosts_col.append(new_host)
hosts_col
#%%

df["hosts"] = hosts_col
df
#%%
_df_save = df.copy()
_df_save.hosts = _df_save.hosts.apply(json.dumps)
_df_save.corr_deseases = _df_save.corr_deseases.apply(json.dumps)
_df_save.to_csv("data/taxons3.csv",index=False)
#%%
# df = pd.read_csv("data/taxons3.csv")
# df.loc[df[df.eppo.isna()].index,"eppo"] = None
# df.loc[df[df.q.isna()].index,"q"] = None
# df.loc[df[df.hosts.isna()].index,"hosts"] = '[]'
# df.loc[df[df.corr_deseases.isna()].index,"corr_deseases"] = '[]'
# df.hosts = df.hosts.apply(json.loads)
# df.corr_deseases = df.corr_deseases.apply(json.loads)
# df

#%%

#@
# По итогу этого всего наш первый пакет данных готов-
# Мы имеем список таксонов патогенов. Каждому в соответствие поставлен: 
# q (гарантировано), 
# eppo code (не гарантировано, если нашёлся в eppo). 
# Список болезней. Которые в свою очередь имеют гарантированный q, т.к взяты из викидаты
# И список хостов. Те в свою очередь имеют гарантированный eppo code  и не гарантированный (если нашёлся) q

#%%

#@

### Второй пакет данных

# На входе имеет таблицу match diseases.
# В ней есть название болезни, патоген, имя хоста (afflict), и eppo code патогена

# Аналогично формируем колонку патогенов. Каждому в соответствие ставятся его eppo код, его болезни, и его аффликты

#%%

df1 = pd.read_csv("match_diseases.csv")
df1
#%%

m_pathogens_u = list(df1["Pathogen"].unique())

_deseases = list([list(df1["Disease"][df1["Pathogen"] == p].unique()) for p in m_pathogens_u])
_eppos = list([list(df1["EPPO code"][df1["Pathogen"] == p].unique())[-1] for p in m_pathogens_u])
_afflicts = list([list(df1["Afflict"][df1["Pathogen"] == p].unique()) for p in m_pathogens_u])

len(m_pathogens_u),len(_deseases),len(_eppos),len(_afflicts)
#%%

# Data from match_diseases in form of taxons dataframe
df2 = pd.DataFrame(data=zip(m_pathogens_u, _eppos, _deseases, _afflicts), columns=["name","eppo","corr_deseases", "afflicts"])
df2
#%%
#@
# Сопоставляем болезни из match_diseases с wikidata

#%%

m_deseases_u:List[str] = []
for deseases in df2["corr_deseases"]:
    m_deseases_u.extend(deseases)
m_deseases_u = np.unique(m_deseases_u)
m_deseases_u
#%%

#@
# У нас уже есть список болезней с wikidata (all_deseases_flat), мы его доставали для первого пакета данных
# Нужно болезни из match_disease сравнить с wikidata

# Это делается при помощи нечёткого поиска библиотеки fuzzywuzzy
# Учитываем метрику схожести >90
# Ставим найденным болезням в соответствие q.


#%%

all_deseases_flat # diseases found in wikidata. tuple (name, q)
#%%

m_desease_name_q_dict:Dict[str,t_q_guaranteed] = {k:v for k,v in all_deseases_flat}
deseases_names = list(m_desease_name_q_dict.keys())
deseases_names

#%%
0/0

@lru_cache(maxsize=None)
def _search_desease_names(name:str, corpus, limit:int):
    res = tuplify(process.extract(name, corpus,limit=limit))
    return res

#%%

def search_desease_name(name:str, corpus:List[str], limit:int=10, threshold=90):
    _res = _search_desease_names(name, tuplify(corpus), limit)
    i = _res[0]
    if i[1]>threshold:
        return i[0]
    else:
        return None
#%%
deseases_search_dict = {d_name: search_desease_name(d_name,deseases_names) for d_name in tqdm(m_deseases_u)}

n_found = sum([bool(i) for i in deseases_search_dict.values()])
print(f'Found {n_found}/{len(m_deseases_u)}')

#%%
match_diseases_dict:Dict[str,Tuple[Optional[str],Optional[t_q]]] = {}


for m_desease_name in m_deseases_u:
    m_desease_name:str
    wd_desease_name:Optional[str] = deseases_search_dict.get(m_desease_name)
    wd_desease_q:Optional[t_q] = m_desease_name_q_dict.get(wd_desease_name)
    match_diseases_dict[m_desease_name] = (wd_desease_name,wd_desease_q) if wd_desease_name else None


#%%
m_deseases_u_dict = match_diseases_dict
#%%
with open("data/disease_corr.json","w+") as fd:
    json.dump(m_deseases_u_dict,fd,indent=2)

#%%

#@
# Найти q для каждого патогена из match_diseases. 
# Это делается поиском по label при помощи sparql запросов, как для хостов в первом пакете

#%%
m_pathogens_wd_search_dict:Dict[str,Optional[t_sparql_response]] = {}

for m_pathogen in tqdm(m_pathogens_u):
    hits = f_q_search_taxon_by_label.cache_info().hits

    res = f_q_search_taxon_by_label(m_pathogen)
    if not res:
        m_pathogens_wd_search_dict[m_pathogen] = None
    else:
        m_pathogens_wd_search_dict[m_pathogen] = res 
    
    if hits == f_q_search_taxon_by_label.cache_info().hits:
        time.sleep(1.5)

len(m_pathogens_wd_search_dict)

#%%
m_pathogens_wd_search_flat_dict = dict()

for k,v in m_pathogens_wd_search_dict.items():
    if not v:
        m_pathogens_wd_search_flat_dict[k] = None
    else:
        flats = [(i["itemLabel"]["value"],i["item"]["value"].split('/')[-1]) for i in v]
        res = None
        if len(flats)>1: # case with multiple findings
            flats_1 = list(filter(lambda x: x[0]==k,flats)) # try find exact match
            if not flats_1:
                res = flats[0]
            else:
                res = flats_1[0]
        else:
            res = flats[0]
                
        m_pathogens_wd_search_flat_dict[k] = res

m_pathogens_wd_search_flat_dict
#%%

df2["q"] = df2.name.map(lambda x: m_pathogens_wd_search_flat_dict.get(x)[1] if m_pathogens_wd_search_flat_dict.get(x) else None)
#%%

df2 = df2[['name','q', 'eppo', 'afflicts', 'corr_deseases']]
df2
#%%
match_diseases_dict:Dict[str,Optional[Tuple[str,t_q]]]
#%%
m_deseases_u_dict_eq = match_diseases_dict
#%%


col_deseases = []
for deseases in df2["corr_deseases"]:
    ds = []
    for desease in deseases:
        if m_deseases_u_dict_eq.get(desease):
            ds.append(m_deseases_u_dict_eq.get(desease))
        else:
            ds.append((desease,None))
    col_deseases.append(ds if ds else None)
col_deseases
#%%
df2["deseases"] = col_deseases
df2

#%%
_df2 = df2[["name","q","eppo","deseases"]]
_df2
#%%
_df_save = df2.copy()
_df_save.afflicts = _df_save.afflicts.apply(json.dumps)
_df_save.corr_deseases = _df_save.corr_deseases.apply(json.dumps)
_df_save.deseases = _df_save.deseases.apply(json.dumps)
_df_save.to_csv("data/taxons6.csv",index=False)
#%%
# df2 = pd.read_csv("data/taxons6.csv")
# df2.loc[df2[df2.eppo.isna()].index,"eppo"] = None
# df2.loc[df2[df2.q.isna()].index,"q"] = None
# df2.loc[df2[df2.afflicts.isna()].index,"afflicts"] = '[]'
# df2.loc[df2[df2.corr_deseases.isna()].index,"corr_deseases"] = '[]'
# df2.loc[df2[df2.deseases.isna()].index,"deseases"] = '[]'
# df2.afflicts = df2.afflicts.apply(json.loads)
# df2.corr_deseases = df2.corr_deseases.apply(json.loads)
# df2.deseases = df2.deseases.apply(json.loads)
# df2
#%%

#%%
#@
# Поиск Afflicts в wikidata

# Unique diseases from match_diseases
m_deseases_u

# Their afflicts dict
m_deseases_u_afflicts = {i:list(df1[df1["Disease"] == i]["Afflict"].unique()) for i in m_deseases_u}

m_afflicts_u = list(df1["Afflict"].unique())

#%%
# Поиск проведён fuzzy_search, результат в c_hosts

#%%
'''

#%%

# Search for afflicts from match_diseases with sparql

hosts_wd_dict = {}
for i in tqdm(m_afflicts_u):
    res = f_q_search_taxon_by_label(i)["results"]["bindings"]
    hosts_wd_dict[i] = res
    time.sleep(1)

#%%

hosts_wd_dict_flat = {}
for k,v in tqdm(hosts_wd_dict.items()):
    if not v: 
        hosts_wd_dict_flat[k] = None
    else:
        hosts_wd_dict_flat[k] = [(i["itemLabel"]["value"], i["item"]["value"].split("/")[-1]) for i in v]

hosts_wd_dict_flat

#%%
hosts_wd_dict_flat_first = {}
for k,v in tqdm(hosts_wd_dict.items()):
    if not v: 
        hosts_wd_dict_flat_first[k] = None
    else:
        hosts_wd_dict_flat_first[k] = (v[0]["itemLabel"]["value"], v[0]["item"]["value"].split("/")[-1])

hosts_wd_dict_flat_first
#%%
hosts_wd_dict_flat_first
#%%
m_afflicts_u
#%%
'''
#%%
# Search for afflicts from match_diseases with API
#%%

# Create the query parameter as a Python dictionary
#%%
0/0

@lru_cache(maxsize=None)
def search_taxon_wikdata(taxon:str):
    url = "https://openrefine-wikidata.toolforge.org/en/api"
    # taxon = "wheat"
    query_param = {
        "query": taxon,
        "limit": 1,
        "properties": [
            {
                "pid": "P31", # instance of
                "v": "Q16521" # taxon
            }
        ]
    }

    # Convert the dictionary to a JSON string
    query_json = json.dumps(query_param)

    # Make the GET request with the query parameter
    response = requests.get(url, params={'query': query_json})

    # Check if the request was successful
    if response.status_code == 200:
        # Parse the JSON response
        data = response.json().get("result")
        if data:
            return (data[0]["name"],data[0]["id"])
        else:
            return None
        
    else:
        print(f"Error: {response.status_code} - {response.text}")
        return None

#%%
# _hosts_wd_api_dict = hosts_wd_api_dict
#%%

hosts_wd_api_dict = {}
for i in tqdm(m_afflicts_u):
    hits = search_taxon_wikdata.cache_info().hits

    res = search_taxon_wikdata(i)
    hosts_wd_api_dict[i] = res
    
    if hits == search_taxon_wikdata.cache_info().hits:
        time.sleep(1)

#%%
'''
a = [[],[],[],[],[]]
for afflict in m_afflicts_u:
    a[0].append(afflict)
    # hosts_wd_dict, hosts_wd_api_dict
    sql_aff = hosts_wd_dict_flat_first[afflict]
    if sql_aff:
        a[1].append(sql_aff[0])
        a[2].append(sql_aff[1])
    else:
        a[1].append(None)
        a[2].append(None)
    api_aff = hosts_wd_api_dict[afflict]
    if api_aff:
        a[3].append(api_aff[0])
        a[4].append(api_aff[1])
    else:
        a[3].append(None)
        a[4].append(None)

df_afflict_wd = pd.DataFrame(zip(*a),columns=["host","sql_name","sql_q","api_name","api_q"])
df_afflict_wd

#%%
df_afflict_wd.to_csv("")
'''

#%%

#@
# По каждому eppo коду патогена извлекаем его хосты. 
# Аналогично пакету 1 находим для них q в Wikidata. 
# Сливаем воедино с Afflict.

#%%

m_eppos_u:List[t_eppo] = list(df1["EPPO code"].unique())

m_eppos_afflict_dict:Dict[t_eppo, List[str]] = {i:list(df1[df1["EPPO code"] == i]["Afflict"].unique()) for i in m_eppos_u}

m_eppos_afflict_dict
#%%

#%%
# По eppo кодам патогенов из match_diseases находим хосты

m_eppo_hosts_dict:Dict[t_eppo_guaranteed,Optional[Dict[Literal['Major host', 'Host'],List[Dict]]]] = dict()
for eppo in tqdm(m_eppos_u):
    hits = get_host.cache_info().hits

    if not eppo:
        m_eppo_hosts_dict[eppo] = None
        continue

    host = get_host(eppo)

    m_eppo_hosts_dict[eppo] = host
    if hits == get_host.cache_info().hits:
        time.sleep(2)

m_eppo_hosts_dict


#%%
m_eppo_hosts_dict_flat:Dict[t_eppo,Optional[List[Tuple[str,t_eppo]]]] = dict()
for eppo,host in m_eppo_hosts_dict.items():

    eppo:t_eppo_guaranteed
    host:Optional[Dict[Literal['Major host', 'Host'],List[Dict]]]    
    
    if not host:
        res = None
    else:
        tuples = []
        for host_item in host.values(): # ['Major host', 'Host']
            for host_i in host_item:
                tuples.append([host_i["full_name"],host_i["eppocode"]])
        res = tuples
    m_eppo_hosts_dict_flat[eppo] = res
    
m_eppo_hosts_dict_flat

#%%
m_flat_hosts_eppo_dict:[Dict[str,t_eppo]] = dict()
for hosts in m_eppo_hosts_dict_flat.values():
    if hosts:
        for host in hosts:
            m_flat_hosts_eppo_dict[host[0]] = host[1]
    
m_flat_hosts_eppo_dict

#%%
# Ищем q хостов

m_flat_hosts_wd_search_dict: Dict[str,Optional[t_sparql_response]] = {}

for m_host in tqdm(m_flat_hosts_eppo_dict):
    m_host:str
    hits = f_q_search_taxon_by_label.cache_info().hits

    res = f_q_search_taxon_by_label(m_host)
    if not res:
        m_flat_hosts_wd_search_dict[m_host] = None
    else:
        m_flat_hosts_wd_search_dict[m_host] = res 
    
    if hits == f_q_search_taxon_by_label.cache_info().hits:
        time.sleep(1.5)

len(m_flat_hosts_wd_search_dict)
#%%

# hosts_names = list(m_flat_hosts_eppo_dict)

# p_size = 30
# l = len(hosts_names)
# pages = int(np.ceil(float(l)/p_size))
# pages
# #%%

# #%%
# m_hosts_wd_search = []

# for i in tqdm(list(range(pages))):
#     items = p(hosts_names,p_size,i)
#     res = f_q_search_taxons_by_labels(items)["results"]["bindings"]
#     m_hosts_wd_search.extend(res)
#     time.sleep(.5)

# len(m_hosts_wd_search)
#%%

m_hosts_wd_search_dict = dict()

for k,v in m_flat_hosts_wd_search_dict.items():
    if not v:
        m_hosts_wd_search_dict[k] = None
    else:
        flats = [(i["itemLabel"]["value"],i["item"]["value"].split('/')[-1]) for i in v]
        res = None
        if len(flats)>1: # case with multiple findings
            flats_1 = list(filter(lambda x: x[0]==k,flats)) # try find exact match
            if not flats_1:
                res = flats[0]
            else:
                res = flats_1[0]
        else:
            res = flats[0]
                
        m_hosts_wd_search_dict[k] = res[1] if res else None

m_hosts_wd_search_dict


#%%
m_flat_hosts_codes_dict:Dict[str,Tuple[Optional[t_q],t_eppo_guaranteed]] = dict()
for host, eppo in m_flat_hosts_eppo_dict.items():
    m_flat_hosts_codes_dict[host] = (m_hosts_wd_search_dict.get(host), eppo)

m_flat_hosts_codes_dict

#%%
m_flat_hosts_codes_triples = [(k,v1,v2) for k,(v1,v2) in m_flat_hosts_codes_dict.items()]

m_flat_hosts_codes_df = pd.DataFrame(m_flat_hosts_codes_triples, columns=["host","q","eppo"])

m_flat_hosts_codes_df.to_csv("data/m_eppo_hosts.csv",index=False)
#%%

m_eppo_hosts_full_dict_flat:Dict[t_eppo,Optional[List[Tuple[str, Optional[t_q],t_eppo_guaranteed]]]] = dict() # like m_eppo_hosts_dict_flat but with qs too

for eppo,hosts in m_eppo_hosts_dict_flat.items():
    new_hosts = []
    if not hosts:
        m_eppo_hosts_full_dict_flat[eppo] = None
        continue
    
    for host in hosts:
        new_hosts.append((host[0],m_flat_hosts_codes_dict[host[0]][0], host[1]))
    
    m_eppo_hosts_full_dict_flat[eppo] = new_hosts

m_eppo_hosts_full_dict_flat

#%%
m_df_afflict_wd = pd.read_csv("c_hosts.csv") # m_name, wd_name, q, eppo
m_df_afflict_wd

#%%

m_afflict_wd_dict:Dict[str,Tuple[str,Optional[t_q],Optional[t_eppo]]] = dict() # m_name: (wd_name, q, eppo)
for i in range(len(m_df_afflict_wd)):
    if str(m_df_afflict_wd.iloc[i,1]) == "nan":
        m_afflict_wd_dict[m_df_afflict_wd.iloc[i,0]] = None
        continue
    m_afflict_wd_dict[m_df_afflict_wd.iloc[i,0]] = (m_df_afflict_wd.iloc[i,1],m_df_afflict_wd.iloc[i,2],m_df_afflict_wd.iloc[i,3])
    
m_afflict_wd_dict
    
#%%
triples:List[Tuple[t_eppo,t_triplet_list,t_triplet_list]] = [] # eppo, own hosts, hosts from afflicts
for eppo in m_eppos_u:
    own_hosts = m_eppo_hosts_full_dict_flat[eppo]
    aff_hosts = [m_afflict_wd_dict[i] for i in m_eppos_afflict_dict[eppo]]
    triples.append((eppo,own_hosts,aff_hosts))
triples
#%%
wd_afflicts:List[t_triplet_list] = []
for afflicts in df2["afflicts"]:
    wds = [m_afflict_wd_dict[afflict] for afflict in afflicts]
    wd_afflicts.append(wds)

df2["wd_afflict"] = wd_afflicts 

df2["eppo_hosts"] = df2["eppo"].map(m_eppo_hosts_full_dict_flat.get)

#%%
# unite eppo_hosts and wd_afflict
hosts = []
for ix,row in df2.iterrows():
    new_host = list(set((row["eppo_hosts"] or []) + (row["wd_afflict"] or [])))
    hosts.append(new_host)

df2["hosts"] = hosts

#%%

#@
# Получаем пакет данных полностью аналогичный первому по структуре. 
# Только в этот раз имеет гарантированные eppo и не гарантированные q

#%%

df2.columns
_df2 = df2.copy()
_df2 = _df2[['name','q', 'eppo', 'hosts', 'deseases'] ]
_df2.columns = ['name','q', 'eppo', 'hosts', 'corr_deseases'] 
_df2
#%%
_df_save = _df2.copy()
_df_save.hosts = _df_save.hosts.apply(json.dumps)
_df_save.corr_deseases = _df_save.corr_deseases.apply(json.dumps)
_df_save.to_csv("data/taxons7.csv",index=False)
#%%
# _df2 = pd.read_csv("data/taxons7.csv")
# _df2.loc[_df2[_df2.eppo.isna()].index,"eppo"] = None
# _df2.loc[_df2[_df2.q.isna()].index,"q"] = None
# _df2.loc[_df2[_df2.hosts.isna()].index,"hosts"] = '[]'
# _df2.loc[_df2[_df2.corr_deseases.isna()].index,"corr_deseases"] = '[]'
# _df2.hosts = _df2.hosts.apply(json.loads)
# _df2.corr_deseases = _df2.corr_deseases.apply(json.loads)
# _df2


#%%
_df = df.copy()
_df = _df[['name','q', 'eppo', 'hosts', 'corr_deseases'] ]

__df = pd.concat([_df2,_df])

#%%
__df.to_csv("a_taxons9.csv", index=False)
#%%
a_df = __df.copy()
#%%

_df_save = a_df.copy()
_df_save.hosts = _df_save.hosts.apply(json.dumps)
_df_save.corr_deseases = _df_save.corr_deseases.apply(json.dumps)
_df_save.to_csv("data/a_taxons9.csv",index=False)
#%%
#@
# Два пакета данных сливаем воедино, a_taxons9
#%%
a_df = pd.read_csv("data/a_taxons9.csv")
a_df.loc[a_df[a_df.eppo.isna()].index,"eppo"] = None
a_df.loc[a_df[a_df.q.isna()].index,"q"] = None
a_df.loc[a_df[a_df.hosts.isna()].index,"hosts"] = '[]'
a_df.loc[a_df[a_df.corr_deseases.isna()].index,"corr_deseases"] = '[]'
a_df.hosts = a_df.hosts.apply(json.loads)
a_df.corr_deseases = a_df.corr_deseases.apply(json.loads)
a_df



#%%

### Третий пакет данных (host_disease_json)

# В первых двух пакетах данных каждому патогена ставятся в соответствия хосты и болезни. 
# Но между хостами и болезнями прямого соответствия нет. 
# Некоторое количество таких соответствий можно получить из match_diseases

#%%
# Table with host-disease correspondance

a_afflicts_u:List[str] = list(df1["Afflict"].unique())
a_afflicts_u
#%%
m_deseases_u
#%%
a_disease_category_dict:Dict[str,List[str]] = {disease:list(df1[df1["Disease"] == disease]['Category'].unique()) for disease in m_deseases_u}
a_disease_category_dict
#%%
a_aflict_diseases_dict:Dict[str,List[Tuple[str,t_q,List[str]]]] = dict()
for afflict in a_afflicts_u:
    _diseases = list(df1[df1["Afflict"] == afflict]['Disease'].unique())
    _new_diseases = []
    
    for _disease in _diseases:
        _category = a_disease_category_dict[_disease]
        if m_deseases_u_dict_eq.get(_disease):
            _t = m_deseases_u_dict_eq[_disease]
            _new_diseases.append((_t[0],_t[1],_category))
        else:
            _new_diseases.append((_disease,None,_category))
    a_aflict_diseases_dict[afflict] = _new_diseases
a_aflict_diseases_dict # host: list of (disease_name, disease_q, list_of_disease_categories)
#%%

p_data:List[Tuple[str, t_q, t_eppo, List[Tuple[str,t_q,List[str]]], t_triplet_list]] = []
for host in a_aflict_diseases_dict:
    _host = m_afflict_wd_dict.get(host) or (host,None,None) # triplet (host, q, eppo)
    _disease = a_aflict_diseases_dict[host] # list of triplets (dis_name, dis_q, dis_cats)
    _pathogens = df1[df1["Afflict"] == host]['Pathogen'].unique()
    _new_pathogens = []  # list of triplets (pathogen, q, eppo)
    for _pathogen in _pathogens:
        _t = m_flat_hosts_codes_dict.get(_pathogen) or (None,None)
        _new_pathogens.append((_pathogen, _t[0], _t[1])) # pathogen, q, eppo
    p_data.append(list(_host)+[list(_disease)]+ [_new_pathogens]) # host, q, eppo, list[(dis_name, dis_q, dis_cats)], pathogens
    
p_data    
#%%
p_df = pd.DataFrame(p_data, columns=["host","q","eppo", "diseases", "pathogens"])
p_df
#%%
#@
# Получаем таблицу такой структуры
# Хосты (аффлиуты)
# Их q и eppo.
# Их список болезней. Для каждой болезни её q(где есть) и категории
# Их список патогенов с q и eppo

#%%
p_df.to_csv("host_disease_1.csv", index=False)
#%%
_df_save = p_df.copy()

_df_save.pathogens = _df_save.pathogens.apply(json.dumps)
_df_save.diseases = _df_save.diseases.apply(json.dumps)
_df_save.to_csv("data/host_disease_1.csv",index=False)

#%%

#######



#%%
#@
### Парсинг таксономии

#%%
# data = pd.read_csv("m_taxons8.csv")
data = pd.read_csv("data/a_taxons9.csv")

#%%
load_f = lambda i: json.loads(repair_json(str(i).replace('(','[').replace(')',']').replace("None","null")) if i else '[]')
#%%

data['hosts']=data['hosts'].apply(load_f)
data['corr_deseases']=data['corr_deseases'].apply(load_f)
data.eppo = data.eppo.apply(lambda x: x or None)
data.q = data.q.apply(lambda x: x or None)
data

#%%
## Uniquify

#%%
_data = data.copy()
df = _data
#%%
len(_data), len(_data["name"].unique()), len(_data.q.unique()), len(_data.eppo.unique())
# (627, 574, 518, 527)

#%%
_df = df
df_by_name = _df.groupby('name').agg({
    'q': 'first',
    'eppo': 'first',
    'hosts': lambda x: list(set(tuplify(i) for sublist in x for i in sublist)),
    'corr_deseases': lambda x: list(set(tuplify(i) for sublist in x for i in sublist))
}).reset_index()

len(_df),len(df_by_name)
#%%
df_by_name
#%%
_df = df_by_name
df_by_q_nulls = _df[_df.q.isnull()]
df_by_q = _df[~_df.q.isnull()].groupby('q').agg({
    'name': 'first',
    'eppo': 'first',
    'hosts': lambda x: list(set(tuple(i) for sublist in x for i in sublist)),
    'corr_deseases': lambda x: list(set(tuple(i) for sublist in x for i in sublist))
}).reset_index()
df_by_q = pd.concat((df_by_q,df_by_q_nulls))
len(_df), len(df_by_q)
#%%
df_by_q
#%%
_df = df_by_q

df_by_eppo_nulls = _df[_df.eppo.isnull()]
df_by_eppo = _df[~_df.eppo.isnull()].groupby('eppo').agg({
    'name': 'first',
    'q': 'first',
    'hosts': lambda x: list(set(tuple(i) for sublist in x for i in sublist)),
    'corr_deseases': lambda x: list(set(tuple(i) for sublist in x for i in sublist))
}).reset_index()
df_by_eppo = pd.concat((df_by_eppo,df_by_eppo_nulls))
len(_df), len(df_by_eppo)
#%%
df_by_eppo

data10 = df_by_eppo.copy()

#%%
_df_save = data10.copy()
_df_save[["name","q","eppo","hosts","corr_deseases"]]
_df_save.hosts = _df_save.hosts.apply(json.dumps)
_df_save.corr_deseases = _df_save.corr_deseases.apply(json.dumps)
_df_save.to_csv("data/a_taxon11.csv",index=False)
#%%

#%%

## End Uniquify

#%%

### Parse pathogens' taxonomy by eppo

#%%

data10 = data10[["name","q","eppo","hosts","corr_deseases"]]
data10

#%%
data12 = data10.copy().reset_index(drop=True)

data12["level"]=[None]* len(data12)
data12["parent"]=[None]* len(data12)
data12 = data12[["level","parent","name","q","eppo","hosts","corr_deseases"]]

for i in tqdm(data12.index):
    
    hits = get_taxonomy.cache_info().hits

    row = data12.loc[i,:]
    eppo = orig_eppo = row["eppo"]
    
    tax_list = get_taxonomy(eppo)
    
    if not tax_list:
        info = get_taxon_info(eppo)
        if info and info.get("replacedby"):
            eppo = info.get("replacedby")
            tax_list = get_taxonomy(eppo)
        else:
            continue
    if not tax_list:
        continue

    for i,t in enumerate(tax_list):
        if i==0:
            parent_ix = None
        else:
            tt = tax_list[i-1]
            parent_eppo = tt["eppocode"]
            parent_ix = data12[data12["eppo"]==parent_eppo].index[0]
        if t["eppocode"] in data12["eppo"].to_list(): # if eppo already in df - just set parent id
            _ix = data12[data12["eppo"]==t["eppocode"]].index
            data12.loc[_ix,"parent"] = parent_ix
            data12.loc[_ix,"level"] = t["level"]
            continue
        else: # else: add new taxon
            l = [t["level"],parent_ix,t["prefname"],None,t["eppocode"],[],[]]
            data12.loc[len(data12)] = l

    if hits == get_taxonomy.cache_info().hits:
        time.sleep(1)
        # pass
#%%
for i in data12.index:
    if "rejected" in data12.loc[i,"name"]:
        _name:str = data12.loc[i,"name"]
        data12.loc[i,"name"] = _name.replace("(rejected name)","").strip()
_df_save = data12.copy()
_df_save.hosts = _df_save.hosts.apply(json.dumps)
_df_save.corr_deseases = _df_save.corr_deseases.apply(json.dumps)
_df_save.to_csv("data/a_taxon12.csv",index=True)
#%%

#%%

## Parse q's for aquired taxonomy

#%%
names = data12[data12["q"].isnull()]["name"]
names

#%%

names_wd_search_dict:Dict[str,Optional[t_sparql_response]] = {}

for name in tqdm(names):
    hits = f_q_search_taxon_by_label.cache_info().hits

    res = f_q_search_taxon_by_label(name)
    if not res:
        names_wd_search_dict[name] = None
    else:
        names_wd_search_dict[name] = res 
    
    if hits == f_q_search_taxon_by_label.cache_info().hits:
        time.sleep(1.5)

len(names_wd_search_dict)

#%%
def _extract_q(x):
    
    ob:list = names_wd_search_dict.get(x)
    if ob is None or len(ob)>1: 
        return None
    
    ob = ob[0]
    q_uri = ob["item"]["value"].split("/")[-1]
    return q_uri
    

extracted_q = data12["name"].apply(_extract_q)
#%%

data13 = data12.copy()

#%%
new_q = [data12.loc[i,"q"] or extracted_q.loc[i] for i in data12["q"].index]
#%%
len(data13[~data13["q"].isnull()]), sum(bool(_i) for _i in new_q)
#%%
data13["q"] = new_q
extracted_q
#%%
_df_save = data13.copy()
_df_save.hosts = _df_save.hosts.apply(json.dumps)
_df_save.corr_deseases = _df_save.corr_deseases.apply(json.dumps)
_df_save.to_csv("data/a_taxon13.csv",index=True)
#%%

## Parse hosts' taxonomy by eppo

#%%
# Gettig flat hosts

hosts_set = set()
hosts_df =pd.DataFrame(columns=["name","q","eppo"])
for hosts in data13["hosts"].to_list():
    for host in hosts:
        host:tuple
        hosts_set.add(host)

for host in hosts_set:
    if not host:
        continue
    if len(host) == 2:
        host = (host[0], None, host[1])
    hosts_df.loc[len(hosts_df)] = host

#%%

# Uniquify hosts

#%%
_df = hosts_df.copy()

len(_df), len(_df["name"].unique()), len(_df.q.unique()), len(_df.eppo.unique())

#%%
df_by_name = _df.groupby('name').agg({
    'q': lambda x: x.mode().min() if x.mode().any() else next((el for el in x if el), None),
    'eppo': lambda x: x.mode().min() if x.mode().any() else next((el for el in x if el), None),
}).reset_index()

len(_df),len(df_by_name)
#%%
df_by_name
#%%
_df = df_by_name

df_by_q_nulls = _df[_df.q.isnull()]
df_by_q = _df[~_df.q.isnull()].groupby('q').agg({
    'name': 'first',
    'eppo': lambda x: x.mode().min() if x.mode().any() else next((el for el in x if el), None),
}).reset_index()
df_by_q = pd.concat((df_by_q,df_by_q_nulls))
len(_df), len(df_by_q)
#%%
df_by_q
#%%

_df = df_by_q

df_by_eppo_nulls = _df[_df.eppo.isnull()]
df_by_eppo = _df[~_df.eppo.isnull()].groupby('eppo').agg({
    'name': 'first',
    'q': lambda x: x.mode().min() if x.mode().any() else next((el for el in x if el), None),
}).reset_index()
df_by_eppo = pd.concat((df_by_eppo,df_by_eppo_nulls))
len(_df), len(df_by_eppo)
#%%
df_by_eppo

hosts_u_df = df_by_eppo[["name","q","eppo",]]
hosts_u_df
#%%

# Get hosts' taxonomy

#%%
hosts_u_df["level"]=[None]* len(hosts_u_df)
hosts_u_df["parent"]=[None]* len(hosts_u_df)
hosts_u_df = hosts_u_df[["level","parent","name","q","eppo"]].reset_index(drop=True)

#%%

for i in tqdm(hosts_u_df.index):
    
    hits = get_taxonomy.cache_info().hits

    row = hosts_u_df.loc[i,:]
    eppo = orig_eppo = row["eppo"]
    
    tax_list = get_taxonomy(eppo)
    
    if not tax_list:
        info = get_taxon_info(eppo)
        if info and info.get("replacedby"):
            eppo = info.get("replacedby")
            tax_list = get_taxonomy(eppo)
        else:
            continue
    if not tax_list:
        continue


    for i,t in enumerate(tax_list):
        if i==0:
            parent_ix = None
        else:
            tt = tax_list[i-1]
            parent_eppo = tt["eppocode"]
            parent_ix = hosts_u_df[hosts_u_df["eppo"]==parent_eppo].index[0]
        if t["eppocode"] in hosts_u_df["eppo"].to_list(): # if eppo already in df - just set parent id
            _ix = hosts_u_df[hosts_u_df["eppo"]==t["eppocode"]].index
            hosts_u_df.loc[_ix,"parent"] = parent_ix
            hosts_u_df.loc[_ix,"level"] = t["level"]
            continue
        else: # else: add new taxon
            l = [t["level"],parent_ix,t["prefname"],None,t["eppocode"]]
            hosts_u_df.loc[len(hosts_u_df)] = l

        
    if hits == get_taxonomy.cache_info().hits:
        time.sleep(.5)
        pass

hosts_u_df
#%%
hosts_u_df.to_csv("plant_taxonomy1.csv")
#%%
# Illiminate duplicated
#%%
vs = hosts_u_df.name.value_counts()
dup_names = vs[vs>1].index.to_list()

for name in dup_names:
    _df = hosts_u_df[hosts_u_df.name==name]
    i = _df.level.argmax()
    
    _s = _df.iloc[i]
    ix = _s.name # index
    res = get_taxon_info(_s.eppo)
    hosts_u_df.loc[ix,"name"] = res["prefname"]

#%%
#%%
hosts_u_df.to_csv("data/plant_taxonomy2_.csv")
#%%

# Parse q's for aquired taxonomy

#%%

names = hosts_u_df[hosts_u_df["q"].isnull()]["name"]
names = [name.strip("x").strip() for name in names]

#%%

hosts_names_wd_search_dict = {}

for name in tqdm(names):
    hits = f_q_search_taxon_by_label.cache_info().hits

    res = f_q_search_taxon_by_label(name)
    if not res:
        hosts_names_wd_search_dict[name] = None
    else:
        hosts_names_wd_search_dict[name] = res 
    
    if hits == f_q_search_taxon_by_label.cache_info().hits:
        time.sleep(1.4)

len(hosts_names_wd_search_dict)
#%%
hosts_names_wd_search_dict
#%%

#%%
def _extract_q(x):
    
    ob:list = hosts_names_wd_search_dict.get(x)
    if ob is None or len(ob)>1: 
        return None
    
    ob = ob[0]
    q_uri = ob["item"]["value"].split("/")[-1]
    return q_uri
    
#%%
extracted_q = hosts_u_df["name"].apply(_extract_q)
#%%
len(hosts_u_df),len(extracted_q)
#%%

names = hosts_u_df[hosts_u_df["q"].isnull()]["name"]
names
#%%
hosts_u_df_2 = hosts_u_df.copy()
hosts_u_df_2
#%%
new_q = [hosts_u_df.loc[i,"q"] or extracted_q.loc[i] for i in hosts_u_df["q"].index]
#%%

hosts_u_df_2["q"] = new_q
#%%
hosts_u_df_2.to_csv("data/plant_taxonomy3.csv")
#%%



#%%

data14 = data13.copy()
vs = data14.name.value_counts()
dup_names = vs[vs>1].index.to_list()
dup_names
#%%
for name in dup_names:
    _df = data14[data14.name==name]
    i = _df.level.argmax()
    
    _s = _df.iloc[i]
    ix = _s.name # index
    res = get_taxon_info(_s.eppo)
    data14.loc[ix,"name"] = res["prefname"]

#%%

_df_save = data14.copy()
_df_save.hosts = _df_save.hosts.apply(json.dumps)
_df_save.corr_deseases = _df_save.corr_deseases.apply(json.dumps)
_df_save.to_csv("data/a_taxon14.csv",index=True)

#%%

vs = data14.name.value_counts()
dup_names = vs[vs>1].index.to_list()
dup_names

#%%

# Link taxons hosts with plant taxonomy

#%%


#%%

taxons_hosts_tuples:List[t_triplet_list] = data14.hosts.apply(tuplify)
taxons_hosts_tuples
#%%
@lru_cache(maxsize=None)
def triplet_compare(t1,t2):
    # _cmp = [t1[0]==t2[0],2*(t1[1]==t2[1]),t1[2]==t2[2]]
    _cmp = [t1[i]==t2[i] for i in range(min(len(t1),len(t2)))]
    # _cmp = (_cmp[0],2*_cmp[1], _cmp[2])
    res = sum(_cmp)+ _cmp[1]
    # return sum((t1[i]==t2[i] for i in range(min(len(t1),len(t2)))))
    # return sum(_cmp)
    return res

#%%

hosts_taxonomy_tuples = [tuplify(r.tolist()) for i,r in hosts_u_df_2[["name","q","eppo"]].iterrows()]
hosts_taxonomy_tuples
#%%

taxons_hosts_tuples_indexes:List[List[int]] = []

ixs = []

for taxon_hosts_list in tqdm(taxons_hosts_tuples):
    taxon_hosts_list:t_triplet_list
    taxons_hosts_indexes = []
    for taxon_triplet in taxon_hosts_list:
        cmp = [triplet_compare(taxon_triplet,host_tuple) for host_tuple in hosts_taxonomy_tuples]
        i = np.argmax(cmp)
        ix = hosts_u_df_2.iloc[i].name # loc index
        if i:
            taxons_hosts_indexes.append(int(ix))
    taxons_hosts_tuples_indexes.append(taxons_hosts_indexes)

taxons_hosts_tuples_indexes
#%%

for i in range(len(taxons_hosts_tuples)):
    if not taxons_hosts_tuples[i]:
        continue 
    
    for ix in taxons_hosts_tuples_indexes[i]:
        # print(hosts_taxonomy_tuples[ix])
        if ix == 0:
            print("Fault",i,ix)
            print(taxons_hosts_tuples[i][ix])
            print(hosts_taxonomy_tuples[ix])
            print()
            print()


#%%
[i for i in taxons_hosts_tuples_indexes if i]
#%%
data15 = data14.copy()
#%%
data15["hosts_ix"] = taxons_hosts_tuples_indexes
#%%
data15 = data15[['level', 'parent', 'name', 'q', 'eppo', 'hosts', 'hosts_ix', 'corr_deseases']]
data15
#%%
data15[data15.hosts.apply(lambda x:bool(x))][["hosts","hosts_ix"]]
#%%

hosts_u_df_2.to_csv("data/plant_taxonomy4.csv")
#%%
data15
#%%

_df_save = data15.copy()
_df_save.hosts = _df_save.hosts.apply(json.dumps)
_df_save.hosts_ix = _df_save.hosts_ix.apply(json.dumps)
_df_save.corr_deseases = _df_save.corr_deseases.apply(json.dumps)
_df_save.to_csv("data/a_taxon15.csv",index=True)

#%%
