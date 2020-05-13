#!/usr/bin/python
import jinja2
import json
import re
import os
import argparse

logger = logging.getLogger()

###########JINJA2 READER#########
TEMPLATE_FILE = "templates/buisciii-infres.j2"


def get_arguments():

        parser = argparse.ArgumentParser(prog = 'common_mash_reference.py', description= 'Search for all mash files and find the representative reference')

        parser.add_argument('-j', '--json', dest="json", metavar="json info file", type=str, required=True, help='REQUIRED. json file with service information.')
        parser.add_argument('-o', '--output', type=str, required=False, help='Output file to save results')

        arguments = parser.parse_args()

        return arguments

    args = get_arguments()

    input_dir = os.path.abspath(args.input_dir)


def main():
    #############FROM API###############################
    #url= 'https://flavia.isciii.es/api/'
    #response = requests.get(url)
    #variables = json.loads(response.text)

    ##################FROM TEXT########################
    templateLoader = jinja2.FileSystemLoader(searchpath="./")
    templateEnv = jinja2.Environment(loader=templateLoader)
    template = templateEnv.get_template(TEMPLATE_FILE)

    with open(args.json, 'r') as service:
        try:
            variables=json.load(service)
        except ValueError:
            print('Decoding JSON has failed')
    #print(variables)
    outputText = template.render(variables)
    file = open(args.output + "/INFRES_"+variables['SERVICE_NUMBER']+".md","wb")
    file.write(outputText.encode("utf-8"))
    file.close()


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        logger.exception(e)
        raise
