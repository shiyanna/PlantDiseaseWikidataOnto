#%%
# Source code for "Automatic Information Retrieval and Extraction Methodology for the Ontology of Plant Diseases"

#%%
import sys
import os
import requests
from SPARQLWrapper import SPARQLWrapper, JSON
import json
import pandas as pd
import numpy as np
from tqdm import tqdm
import time
from functools import lru_cache
from bs4 import BeautifulSoup

#%%
endpoint_url = "https://query.wikidata.org/sparql"
#%%
#=========
# EPPO


#%%

# EPPO api token. Get it here: https://data.eppo.int
_AUTHTOKEN = 'your-token-here'

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

# Search list of taxons in eppo
def search_taxons(taxons:list,authtoken = None):
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

#%%
0/0 # Be careful not to reexecute to not lose cache

# Get taxon's hosts from eppo by eppo code
@lru_cache(maxsize=None)
def get_host(taxon,authtoken = None):
    if not authtoken:
        global _AUTHTOKEN
        authtoken = _AUTHTOKEN
        
    params = {
        'authtoken': authtoken
    }

    url = 'https://data.eppo.int/api/rest/1.0/taxon/%s/hosts'
    
    res = json_get_request(url=url%taxon, params=params)
    return res

# END EPPO
#=========

#%%

#%%
0/0 # Be careful not to reexecute to not lose cache

# Execute sparql query in wikidata
@lru_cache(maxsize=None)
def get_results(query, endpoint_url = None):
    if not endpoint_url:
        endpoint_url = "https://query.wikidata.org/sparql"
    user_agent = "WDQS-example Python/%s.%s" % (sys.version_info[0], sys.version_info[1])
    sparql = SPARQLWrapper(endpoint_url, agent=user_agent)
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    return sparql.query().convert()

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
subdeseases_qs = [i[1] for i in subdeseases_tuples]
subdeseases_qs.append("Q2662845")
subdeseases_qs
#%%

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

res = f_q_all_deseases(subdeseases_qs)
all_deseases = res["results"]["bindings"]
len(all_deseases)

#%%

os.makedirs("data", exist_ok=True)

#%%
with open("data/all_deseases.json", "w+") as fd:
    json.dump(all_deseases, fd, indent=2)
#%%
all_deseases_flat = [(i["itemLabel"]["value"],i["item"]["value"].split("/")[-1]) for i in all_deseases]
all_deseases_flat
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
causes

#%%

with open("data/all_deseases.json", "w+") as fd:
    json.dump(causes, fd, indent=2)
#%%

causes_flat = []

for item in causes["results"]["bindings"]:
    q = item["caused"]["value"].split("/")[-1]
    label = item["causedLabel"]["value"]
    causes_flat.append((q,label))

deseases_causes_flat = []

for item in causes["results"]["bindings"]:
    desease_label = item["itemLabel"]["value"]
    desease_q = item["item"]["value"].split("/")[-1]
    q = item["caused"]["value"].split("/")[-1]
    label = item["causedLabel"]["value"]
    deseases_causes_flat.append((desease_q,desease_label,q,label))

#%%
df_deseases = pd.DataFrame(data=deseases_causes_flat, columns=["desease_q","desease_label","cause_q","cause_label"])
df_deseases
#%%
#%%
causes_dict = {label:q for q,label in causes_flat}
names = list(causes_dict)

taxons = search_taxons(names)

causes_triples = [[names[i],causes_dict[names[i]], taxons[i]] for i in range(len(causes_dict))]
causes_triples

df = pd.DataFrame(data = causes_triples, columns=["name", "q", "eppo"])
df

#%%

# Get hosts from eppo. Makes api call, so can be interrupted on rate limit

hosts = []
for eppo in tqdm(df["eppo"]):
    # print(eppo)
    if not eppo:
        hosts.append(None)
        continue
    host = get_host(eppo)
    # print(host)
    hosts.append(host)
    time.sleep(2)

hosts
#%%

hosts_tuples = []

for host in hosts:
    if not host:
        res = None
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
hosts = df['hosts']

hosts_flat = []
for host in hosts:
    if host:
        hosts_flat.extend(host)

#%%
hosts_dict = {k:v for k,v in hosts_flat}

hosts_tuples = list(hosts_dict.items())

hosts_dict_names = {v:k for k,v in hosts_dict.items()}

hosts_names = [i[1] for i in hosts_tuples]

#%%
# df.to_csv("taxons.csv",index=False)
#%%
corr_deseases = []
for q in tqdm(df["q"]):
    _df = df_deseases[df_deseases["cause_q"]==q][["desease_q","desease_label"]]
    corr_desease = [(a,b) for a,b in zip(_df["desease_label"],_df["desease_q"])]
    corr_deseases.append(corr_desease)
corr_deseases

df["corr_deseases"] = corr_deseases
df
#%%

q_search_by_label = '''
SELECT distinct ?item ?itemLabel WHERE{  
  ?item ?label ?LAB .  
  ?article schema:about ?item .
  VALUES ?LAB { %s }
  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }    
}
'''

def f_q_search_by_labels(labels): 
    if type(labels) == list:
        s_labels = " ".join([f'"{label}"@en' for label in labels])
    else:
        s_labels = f'"{labels}"@en'
    # print(q_search_by_label%s_labels)
    res = get_results(q_search_by_label%s_labels)
    return res

#%%

def p(items, size, page): # paginator
    return items[page*size:(page+1)*size]
#%%

p_size = 30
l = len(hosts_names)
pages = int(np.ceil(float(l)/p_size))
pages

#%%
wd_items = []
for i in tqdm(list(range(pages))):
    items = p(hosts_names,p_size,i)
    res = f_q_search_by_labels(items)["results"]["bindings"]
    wd_items.extend(res)
    time.sleep(.5)

len(wd_items)

#%%

wd_items_dict = {item["itemLabel"]["value"]:item["item"]["value"].split('/')[-1] for item in wd_items}
wd_items_dict

#%%

hosts_dict_complete = {name:[hosts_dict_names[name], wd_items_dict.get(name)] for name in hosts_dict_names}
hosts_dict_complete

#%%

# Appending q found in wikidata
new_hosts = []
for host in df["hosts"]:
    if not host:
        new_hosts.append(None)
        continue
    new_host = []
    for host_i in host:
        name = host_i[1]
        eppo = host_i[0]
        q = hosts_dict_complete[name][1]
        triple = (name,q,eppo)
        new_host.append(triple)
    new_hosts.append(new_host)
new_hosts
#%%

df["hosts"] = new_hosts
df
#%%
# Formating
#%%
cols = ["hosts","corr_deseases"]
_df = df.copy()

for col in cols:
    new_col = df[col].map(lambda x:json.dumps(x) if x else None)
    new_col
    _df[col] = new_col
_df
#%%

_df.to_csv("taxons3.csv",index=False)


#%%
@lru_cache(maxsize=None)
def get_datasheet(eppo_code):
    if not eppo_code: 
        return None
    url = f"https://gd.eppo.int/taxon/{eppo_code}/datasheet"
    response = requests.get(url)
    soup = BeautifulSoup(response.text, 'html.parser')

    datasheet_content = soup.find("div", class_="datasheet")
    return datasheet_content.text if datasheet_content else None

#%%

def get_raw_datasheet(eppo_code):
    if not eppo_code: 
        return None
    url = f"https://gd.eppo.int/taxon/{eppo_code}/datasheet"
    response = requests.get(url)
    
    text = response.text
    text = text.replace("display: none;","display: true;")
    
    return text

#%%

soup = BeautifulSoup(response.text, 'html.parser')
soup

#%%
datasheet_content = soup.find("div", class_="datasheet")
datasheet_content

#%%

datasheets = []

for eppo in tqdm(df["eppo"]):
    hits = get_datasheet.cache_info().hits
    datasheet = get_datasheet(eppo)
    datasheets.append(datasheet)
    if hits == get_datasheet.cache_info().hits:
        time.sleep(2)

#%%
s1 =  "DETECTION AND IDENTIFICATION"
s2 =  "Symptoms"

valid_datasheets = [(datasheet if datasheet and s1 in datasheet else None) for datasheet in datasheets]

valid_datasheets1 = [(datasheet if datasheet and s2 in datasheet else None) for datasheet in valid_datasheets]

df["datasheet"] = valid_datasheets1
#%%

_df = df.copy()

for col in cols:
    new_col = df[col].map(lambda x:json.dumps(x) if x else None)
    new_col
    _df[col] = new_col
_df
#%%

_df.to_csv("taxons4.csv",index=False)

#%%

df1 = pd.read_csv("match_diseases.csv")
df1
#%%

m_pathogens_u = list(df1["Pathogen"].unique())

_deseases = list([list(df1["Disease"][df1["Pathogen"] == p].unique()) for p in m_pathogens_u])

_eppos = list([list(df1["EPPO code"][df1["Pathogen"] == p].unique())[-1] for p in m_pathogens_u])

_afflicts = list([list(df1["Afflict"][df1["Pathogen"] == p].unique()) for p in m_pathogens_u])

#%%

# Data from match_deseases in form of taxons dataframe
df2 = pd.DataFrame(data=zip(m_pathogens_u, _eppos, _deseases, _afflicts), columns=["name","eppo","corr_deseases", "afflicts"])
df2
#%%

m_deseases_u = []
for deseases in df2["corr_deseases"]:
    m_deseases_u.extend(deseases)
m_deseases_u = np.unique(m_deseases_u)
m_deseases_u
#%%

all_deseases_flat # diseases found in wikidata. tuple (name, q)
#%%

def str_dist(s1:str,s2:str): # words intersection over union
    stuff = ["'s","(",")","`","ʻ","'",'"','.']
    
    ss1 = s1.replace("-"," ").lower()
    for s in stuff:
        ss1 = ss1.replace(s,"")
    set1 = set(ss1.split())
    
    ss2 = s2.replace("-"," ").lower()
    for s in stuff:
        ss2 = ss2.replace(s,"")
    set2 = set(ss2.split())
    
    jaccard_index = len(set1.intersection(set2)) / len(set1.union(set2))
    return jaccard_index

#%%
# Compare each desease from match_deseases.csv (m_deseases_u) to deseases from wikidata (all_deseases_flat)
# Structure:  m_disease_name: (score, wd_disease_name, q)

threshold = .5

m_deseases_u_dict = dict()
for desease in m_deseases_u:
    indexes = [(round(str_dist(i[0],desease),2),i[0],i[1]) for i in all_deseases_flat if str_dist(i[0],desease)>=threshold]
    if indexes:
        m_deseases_u_dict[str(desease)] = sorted(indexes,key=lambda x:x[0],reverse=True)

#%%
print('Total unique diseases in "match_diseases.csv":', len(m_deseases_u))
print(f'Found matches above {threshold} threshold',len(m_deseases_u_dict))
#%%
m_deseases_u_dict
#%%
with open("disease_corr.json","w+") as fd:
    json.dump(m_deseases_u_dict,fd,indent=2)

#%%
#%%
# Find perfect matches

m_deseases_u_dict_eq = dict()

for k,v in m_deseases_u_dict.items():
    if v[0][0]==1.0:
        m_deseases_u_dict_eq[k] = (v[0][1],v[0][2])

m_deseases_u_dict_eq
#%%
print('Perfect matched diseases from "match_diseases.csv to wikidata:', len(m_deseases_u_dict_eq))
#%%
# Search for pathogens from match_diseases

#%%

p_size = 30
l = len(m_pathogens_u)
pages = int(np.ceil(float(l)/p_size))
pages
#%%

q_search_taxon_by_label = '''
SELECT distinct ?item ?itemLabel WHERE{  
  ?item ?label ?LAB .  
  ?item wdt:P31 wd:Q16521 . # # p31 instance of  q16521 taxon
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

@lru_cache(maxsize=None)
def f_q_search_taxon_by_label(label:str): 
    s_labels = f'"{label}"@en'
    res = get_results(q_search_taxon_by_label%s_labels)
    return res

#%%
m_pathogens_wd_search_dict = {}

for m_pathogen in tqdm(m_pathogens_u):
    hits = f_q_search_taxon_by_label.cache_info().hits

    res = f_q_search_taxon_by_label(m_pathogen)["results"]["bindings"]
    if not res:
        m_pathogens_wd_search_dict[m_pathogen] = None
    else:
        m_pathogens_wd_search_dict[m_pathogen] = res 
    
    if hits == f_q_search_taxon_by_label.cache_info().hits:
        time.sleep(2)

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

#%%

col_deseases = []
for deseases in df2["corr_deseases"]:
    ds = []
    for desease in deseases:
        if desease in m_deseases_u_dict_eq:
            ds.append(m_deseases_u_dict_eq.get(desease))
        else:
            ds.append((desease,None))
    col_deseases.append(ds if ds else None)
col_deseases

df2["deseases"] = col_deseases
df2

#%%
_df2 = df2[["name","q","eppo","deseases"]]
_df2
#%%
_df2.to_csv("taxons6.csv")
#%%
df2.columns
#%%

# Unique diseases from match_diseases
# m_deseases_u = list(df1["Disease"].unique())
m_deseases_u

# Their afflicts dict
m_deseases_u_afflicts = {i:list(df1[df1["Disease"] == i]["Afflict"].unique()) for i in m_deseases_u}

m_afflicts_u = list(df1["Afflict"].unique())

#%%
# This was done once and then fixed manualy
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
# Exp
# Afflict for eppo, for comparing afflict with eppo host
#%%

m_eppos_u = list(df1["EPPO code"].unique())

m_eppos_afflict_dict = {i:list(df1[df1["EPPO code"] == i]["Afflict"].unique()) for i in m_eppos_u}

m_eppos_afflict_dict
#%%


m_afflicts_u = list(df1["Afflict"].unique())

m_afflict_eppos_dict = {i:list(df1[df1["Afflict"] == i]["EPPO code"].unique()) for i in m_afflicts_u}

m_afflict_eppos_dict

#%%
len(m_eppos_u), len(m_eppos_afflict_dict)
#%%

m_eppo_hosts_dict = dict()
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
m_eppo_hosts_dict_flat = dict()
for eppo,host in m_eppo_hosts_dict.items():
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
m_flat_hosts_eppo_dict = dict()
for hosts in m_eppo_hosts_dict_flat.values():
    if hosts:
        for host in hosts:
            m_flat_hosts_eppo_dict[host[0]] = host[1]
    
m_flat_hosts_eppo_dict

#%%

m_flat_hosts_wd_search_dict = {}

for m_host in tqdm(m_flat_hosts_eppo_dict):
    hits = f_q_search_taxon_by_label.cache_info().hits

    res = f_q_search_taxon_by_label(m_host)["results"]["bindings"]
    if not res:
        m_flat_hosts_wd_search_dict[m_host] = None
    else:
        m_flat_hosts_wd_search_dict[m_host] = res 
    
    if hits == f_q_search_taxon_by_label.cache_info().hits:
        time.sleep(2)

len(m_flat_hosts_wd_search_dict)
#%%

hosts_names = list(m_flat_hosts_eppo_dict)

p_size = 30
l = len(hosts_names)
pages = int(np.ceil(float(l)/p_size))
pages
#%%

#%%
m_hosts_wd_search = []

for i in tqdm(list(range(pages))):
    items = p(hosts_names,p_size,i)
    res = f_q_search_taxons_by_labels(items)["results"]["bindings"]
    m_hosts_wd_search.extend(res)
    time.sleep(.5)

len(m_hosts_wd_search)
#%%

m_hosts_wd_search_dict = {item["itemLabel"]["value"]:item["item"]["value"].split('/')[-1] for item in m_hosts_wd_search}
m_hosts_wd_search_dict
#%%
m_flat_hosts_codes_dict = dict()
for host, eppo in m_flat_hosts_eppo_dict.items():
    m_flat_hosts_codes_dict[host] = (m_hosts_wd_search_dict.get(host), eppo)

m_flat_hosts_codes_dict

#%%
m_flat_hosts_codes_triples = [(k,v1,v2) for k,(v1,v2) in m_flat_hosts_codes_dict.items()]

m_flat_hosts_codes_df = pd.DataFrame(m_flat_hosts_codes_triples, columns=["host","q","eppo"])

m_flat_hosts_codes_df.to_csv("m_all_hosts.csv",index=False)
#%%

m_eppo_hosts_full_dict_flat = dict() # like m_eppo_hosts_dict_flat but with qs too

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

m_afflict_wd_dict = dict() # m_name: (wd_name, q, eppo)
for i in range(len(m_df_afflict_wd)):
    if str(m_df_afflict_wd.iloc[i,1]) == "nan":
        m_afflict_wd_dict[m_df_afflict_wd.iloc[i,0]] = None
        continue
    m_afflict_wd_dict[m_df_afflict_wd.iloc[i,0]] = (m_df_afflict_wd.iloc[i,1],m_df_afflict_wd.iloc[i,2],m_df_afflict_wd.iloc[i,3])
    
m_afflict_wd_dict
    
#%%
triples = [] # eppo, own hosts, hosts from afflicts
for eppo in m_eppos_u:
    own_hosts = m_eppo_hosts_full_dict_flat[eppo]
    aff_hosts = [m_afflict_wd_dict[i] for i in m_eppos_afflict_dict[eppo]]
    triples.append((eppo,own_hosts,aff_hosts))

#%%
wd_afflicts = []
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
    # print(row)

df2["hosts"] = hosts
#%%

df2.columns
_df2 = df2.copy()
_df2 = _df2[['name','q', 'eppo', 'hosts', 'deseases'] ]
_df2.columns = ['name','q', 'eppo', 'hosts', 'corr_deseases'] 

cols = ["hosts","corr_deseases"]

for col in cols:
    new_col = _df2[col].map(lambda x:json.dumps(x) if x else None)
    new_col
    _df2[col] = new_col

# _df2.to_csv("m_taxons8.csv", index=False)
#%%
_df = df.copy()
_df = _df[['name','q', 'eppo', 'hosts', 'corr_deseases'] ]

for col in cols:
    new_col = _df[col].map(lambda x:json.dumps(x) if x else None)
    new_col
    _df[col] = new_col
#%%
#%%
__df = pd.concat([_df2,_df])

#%%
__df.to_csv("a_taxons9.csv", index=False)

a_df = __df
#%%


#%%
# Table with host-disease correspondance

a_afflicts_u = list(df1["Afflict"].unique())
a_afflicts_u
#%%
m_deseases_u
#%%
a_disease_category_dict = {disease:list(df1[df1["Disease"] == disease]['Category'].unique()) for disease in m_deseases_u}
a_disease_category_dict
#%%
a_aflict_diseases_dict = dict()
for afflict in a_afflicts_u:
    _diseases = list(df1[df1["Afflict"] == afflict]['Disease'].unique())
    _new_diseases = []
    
    for _disease in _diseases:
        _category = a_disease_category_dict[_disease]
        if _disease in m_deseases_u_dict_eq:
            _t = m_deseases_u_dict_eq[_disease]
            _new_diseases.append((_t[0],_t[1],_category))
        else:
            _new_diseases.append((_disease,None,_category))
    a_aflict_diseases_dict[afflict] = _new_diseases
a_aflict_diseases_dict # host: list of (disease_name, disease_q, list_of_disease_categories)
#%%

p_data = []
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
p_df.to_csv("host_disease_1.csv", index=False)
#%%
cols = ["diseases","pathogens"]

_df = p_df.copy()
for col in cols:
    new_col = p_df[col].map(lambda x:json.dumps(x) if x else None)
    new_col
    _df[col] = new_col
_df

_df.to_csv("host_disease_json_1.csv", index=False)

#%%

#######


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

#%%
# data = pd.read_csv("m_taxons8.csv")
data = pd.read_csv("a_taxons9.csv")
#%%

#%%
load_f = lambda i: json.loads(repair_json(str(i).replace('(','[').replace(')',']').replace("None","null")) if i else '[]')
#%%

data['hosts']=data['hosts'].apply(load_f)
data['corr_deseases']=data['corr_deseases'].apply(load_f)
data.eppo = data.eppo.apply(lambda x: x or None)
data.q = data.q.apply(lambda x: x or None)

#%%
## Populate eppos

@lru_cache(maxsize=None)
def read_local_eppo():
    _df1 = pd.read_csv("eppo/gafname.txt")
    _df2 = pd.read_csv("eppo/gainame.txt")
    _df3 = pd.read_csv("eppo/pflname.txt")
    _df = pd.concat((_df1,_df2,_df3))    
    _names = _df["fullname"].dropna().tolist()
    return _df, _names
#%%
#%%
@lru_cache(maxsize=None)
def search_local_eppo_names(name, limit = 20):
    if name in _dfindings:
        return _dfindings[name]
    _df,_names = read_local_eppo()
    column_as_list = _names
    res = tuplify(process.extract(name, column_as_list,limit=limit))
    return res

#%%

df_missing_eppo = data[data.eppo.isnull()]
misses = df_missing_eppo["name"]
len(misses)

#%%

findings = []

for st in tqdm(misses):
    f = search_local_eppo_names(st)
    findings.append((st,f))

len(findings)
#%%
findings
#%%
found = dict()
not_found = set()
for f in findings:
    trg = f[1][0]
    if trg[1]>90:
        found[f[0]] = trg[0]
    else:
        not_found.add(f[0])
found, not_found
#%%
len(found), len(not_found)

#%%
_eppo_df,_names = read_local_eppo()
#%%

print(len(data[data.eppo.isnull()]))

for k,v in found.items():
    ix = data[data.name==k].index
    new_eppo = _eppo_df[_eppo_df.fullname==v].code.iloc[0]
    data.loc[ix,"eppo"] = new_eppo
    
print(len(data[data.eppo.isnull()]))


#%%
## Uniquify

#%%
_data = data.copy()
df = _data
#%%
len(_data), len(_data["name"].unique()), len(_data.q.unique()), len(_data.eppo.unique())
#%%
tuplify = lambda arr: tuple(i if type(i) == str else tuple(i) for i in arr if i) if arr else ()
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
df_by_eppo.to_csv('a_taxon11.csv',index=False)
#%%

## End Uniquify
#%%
_df_save = data10.copy()
_df_save.columns
#%%
_df_save.hosts = _df_save.hosts.apply(json.dumps)
_df_save.corr_deseases = _df_save.corr_deseases.apply(json.dumps)
#%%
_df_save = _df_save[["name","q","eppo","hosts","corr_deseases"]]
_df_save
#%%
_df_save.to_csv('a_taxon11.csv',index=False)
#%%

## Parse pathogens' taxonomy by eppo

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


data12.to_csv("a_taxon12.csv")

#%%

## Parse q's for aquired taxonomy

#%%
names = data12[data12["q"].isnull()]["name"]
names

#%%

names_wd_search_dict = {}

for name in tqdm(names):
    hits = f_q_search_taxon_by_label.cache_info().hits

    res = f_q_search_taxon_by_label(name)["results"]["bindings"]
    if not res:
        names_wd_search_dict[name] = None
    else:
        names_wd_search_dict[name] = res 
    
    if hits == f_q_search_taxon_by_label.cache_info().hits:
        time.sleep(2)

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

## Parse hosts' taxonomy by eppo

#%%
# Gettig flat hosts

hosts_set = set()
hosts_df =pd.DataFrame(columns=["name","q","eppo"])
for hosts in data12["hosts"].to_list():
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

# Parse q's for aquired taxonomy

#%%

names = hosts_u_df[hosts_u_df["q"].isnull()]["name"]
names = [name.strip("x").strip() for name in names]

#%%

hosts_names_wd_search_dict = {}

for name in tqdm(names):
    hits = f_q_search_taxon_by_label.cache_info().hits

    res = f_q_search_taxon_by_label(name)["results"]["bindings"]
    if not res:
        hosts_names_wd_search_dict[name] = None
    else:
        hosts_names_wd_search_dict[name] = res 
    
    if hits == f_q_search_taxon_by_label.cache_info().hits:
        time.sleep(2)

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
hosts_u_df_2.to_csv("plant_taxonomy2.csv")
#%%



#%%



#%%



#%%



#%%

