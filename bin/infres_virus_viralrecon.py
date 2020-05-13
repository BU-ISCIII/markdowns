#!/usr/bin/python
import jinja2
import json
import re
import os

###########JINJA2 READER#########
templateLoader = jinja2.FileSystemLoader(searchpath="./")
templateEnv = jinja2.Environment(loader=templateLoader)
TEMPLATE_FILE = "template/INFRES_VIRUS_VIRALRECON.j2"
template = templateEnv.get_template(TEMPLATE_FILE)

#############FROM API###############################
#url= 'https://flavia.isciii.es/api/'
#response = requests.get(url)
#variables = json.loads(response.text)

##################FROM TEXT########################
with open("data/service_data_mapping.json", 'r') as service:
    try:
        variables=json.load(service)
    except ValueError:
        print('Decoding JSON has failed')
#print(variables)
outputText = template.render(variables)
file = open("results/INFRES_"+variables['SERVICE_NUMBER']+".md","wb")
file.write(outputText.encode("utf-8"))
file.close()
