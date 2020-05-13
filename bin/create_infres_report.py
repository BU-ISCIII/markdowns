#!/usr/bin/python
import jinja2
import json
import re
import os
import argparse
import datetime

logger = logging.getLogger()

###########JINJA2 READER#########
TEMPLATE_FILE = "templates/buisciii-infres.j2"


def get_arguments():

        parser = argparse.ArgumentParser(prog = 'common_mash_reference.py', description= 'Search for all mash files and find the representative reference')

        parser.add_argument('-j', '--json', dest="json", metavar="json info file", type=str, required=True, help='REQUIRED. json file with service information.')
        parser.add_argument('-o', '--output', type=str, required=False, help='Output file to save results')

        arguments = parser.parse_args()

        return arguments


def main():
    # get arg
    args = get_arguments()

    #Create log file with date and time
    right_now = str(datetime.date.today())
    right_now_full = "_".join(right_now.split(" "))

    log_filename = 'infres_creation' + "_" + right_now_full + ".log"
    log_full_path = os.path.join(args.output, log_filename)

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s:%(message)s')

    file_handler = logging.FileHandler(log_full_path)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    #stream_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)


    #####################START PIPELINE################

    logger.info(args)

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
        except ValueError as e:
            logger.error("json decoding has failed")
            logger.exception(e)

    #print(variables)
    outputText = template.render(variables)
    file = open(args.output + "/INFRES_"+variables['SERVICE_NUMBER']+".md","wb")
    file.write(outputText.encode("utf-8"))
    file.close()

    logger.info("Infres has being successfully rendered.")

if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        logger.exception(e)
        raise
