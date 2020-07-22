#!/usr/bin/python
import jinja2
import json
import re
import os
import argparse
import datetime
import logging
import pymysql
import sys

logger = logging.getLogger()

###########JINJA2 READER#########
BASEPATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_FILE = "templates/buisciii-infres.j2"


def get_arguments():
    parser = argparse.ArgumentParser(prog = 'create_infres_report.py', description= 'Create markdown file of services for report')
    parser.add_argument('-j', '--json', dest="json", metavar="json info file", type=str, required=False, help='REQUIRED if the request to the database is not being done. json file with service information.')
    parser.add_argument('-o', '--output', type=str, required=False, help='Output markdown file to save results')
    parser.add_argument('-s', '--service', type=str, required=False, help='Service ID to make the request to the database. REQUIRED if request is being done.')
    parser.add_argument('-m', '--samples', type=str, required=False, help='Path to the samples ID file to filter the samples in the request. REQUIRED if request to databse is being done.')
    parser.add_argument('-f', '--save_json', type=str, required=False, help='Path to save the output json file. REQUIRED if request to databse is being done.')

    arguments = parser.parse_args()

    return arguments

def check_samples(sample,samples,sample_list):
    sample_test_rename = sample.split("_",1)[0] #We split the sample_id in the project by "_" because we rename them.
    if sample in sample_list or sample_test_rename in sample_list: #We check if the sample is in the samples_id list
        if sample not in samples: #We check if the samples is already in the sample list
            samples.append(sample) #If it is in samples_id and it's not int sample list, we add it
    return samples


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
    logger.info("Starting infres report creation...")
    logger.info(args)
    if args.samples:
        with open(args.samples) as sample_file:
            sample_list = [line.rstrip() for line in sample_file]

    ##################FROM TEXT########################
    #print(BASEPATH)
    templateLoader = jinja2.FileSystemLoader(searchpath=BASEPATH)
    templateEnv = jinja2.Environment(loader=templateLoader)
    template = templateEnv.get_template(TEMPLATE_FILE)

    if args.json:
        with open(args.json, 'r') as service:
            try:
                json_data=json.load(service)
            except ValueError as e:
                logger.error("json decoding has failed")
                logger.exception(e)
    else:
        if not args.service:
            logger.error("When json file not provided, you must provide a service id. --service")
            exit()

        #############FROM API###############################
        # Open database connection
        db = pymysql.connect("flavia.isciii.es","django","django77","iSkyLIMS" )

        # prepare a cursor object using cursor() method
        cursor = db.cursor()

        resolution_query="SELECT resolution.resolutionNumber 'Resolution Number' FROM iSkyLIMS_drylab_service service, iSkyLIMS_drylab_resolution resolution WHERE resolution.resolutionServiceID_id = service.id and service.serviceRequestNumber='"+args.service+"' ORDER BY resolution.resolutionNumber DESC LIMIT 1;"
        cursor.execute(resolution_query)
        resolution_data=cursor.fetchall()
        last_resolution=resolution_data[0][0]


        # execute SQL query using execute() method.
        query="SELECT service.serviceRequestNumber 'Service Number', resolution.resolutionNumber 'Resolution Number', resolution.resolutionFullNumber 'Service Name', service.serviceSeqCenter 'Sequencing Center', service.serviceRunSpecs 'Run Specifications', service.serviceCreatedOnDate 'Request Date', resolution.resolutionDate 'Resolution Date', resolution.resolutionEstimatedDate 'Estimated Delivery Date', delivery.deliveryDate 'Delivery Date', resolution.resolutionOnQueuedDate 'Queue Date', resolution.resolutionOnInProgressDate 'In progress Date', service.serviceNotes 'Notes', users.username 'Username', users.first_name 'First Name', users.last_name 'Last Name', profile.profilePosition 'Position', center.centerName 'Center', profile.profileArea 'Area', profile.profileExtension 'Phone Extension', users.email 'E-mail', sample.sampleName 'SampleName', project.projectName 'Project Name', runname.runName 'Run Name' FROM iSkyLIMS_drylab_service service, iSkyLIMS_drylab_service_serviceProjectNames bridge, iSkyLIMS_wetlab_projects project, iSkyLIMS_wetlab_runprocess runname, iSkyLIMS_wetlab_samplesinproject sample, iSkyLIMS_drylab_resolution resolution, iSkyLIMS_drylab_delivery delivery, auth_user users, django_utils_profile profile, django_utils_center center WHERE bridge.service_id = service.id and bridge.projects_id = project.id and sample.project_id_id = project.id and runname.id = project.runprocess_id_id and resolution.resolutionServiceID_id = service.id and delivery.deliveryResolutionID_id = resolution.id and service.serviceUserId_id = users.id and users.id = profile.profileUserID_id and profile.profileCenter_id = center.id and resolution.resolutionNumber='"+last_resolution+"';"

        cursor.execute(query)
        data = cursor.fetchall()

        # Fetch a single row using fetchone() method.s
        json_data = {}
        service={}
        user={}
        runs=[]

        for row in data:
            if 'SERVICE_NUMBER' not in json_data:
                json_data["SERVICE_NUMBER"] = row[0]
                json_data["RESOLUTION_NUMBER"] = row[1]
                service["ID"] = row[2]
                service["SEQUENCING_CENTER"] = row[3]
                service["RUN_SPECIFICATIONS"] = row[4]
                req_date=str(row[5].year)+'/'+str(row[5].month)+'/'+str(row[5].day)
                service["REQUEST_DATE"] = req_date
                resol_date=str(row[6].year)+'/'+str(row[6].month)+'/'+str(row[6].day)
                service["RESOLUTION_DATE"] = resol_date
                prog_date=str(row[10].year)+'/'+str(row[10].month)+'/'+str(row[10].day)
                service["IN_PROGRESS_DATE"] = prog_date
                est_date=str(row[7].year)+'/'+str(row[7].month)+'/'+str(row[7].day)
                service["ESTIMATED_DELIVERY_DATE"] = est_date
                deliv_date=str(row[8].year)+'/'+str(row[8].month)+'/'+str(row[8].day)
                service["DELIVERY_DATE"] = deliv_date
                notes = row[11]
                notes = notes.replace('\r', '')
                notes = notes.replace('\n', ' ')
                service["SERVICE_NOTES"] = notes
                json_data["SERVICE"] = service
                user["USERNAME"] = row[12]
                user["FIRST_NAME"] = row[13]
                user["LAST_NAME"] = row[14]
                user["POSITION"] = row[15]
                user["CENTER"] = row[16]
                user["UNIT"] = row[17]
                if row[18] != "-":
                    user["PHONE"] = row[18]
                user["EMAIL"] = row[19]
                json_data["USER"] = user
            test_run=next((run for run in runs if run["RUN_NAME"] == row[22]), False) #Check if the run exist
            if test_run != False: #if run already exist
                project={}
                test_poject=next((project for project in projects if project["PROJECT_NAME"] == row[21]), False) #Check if run exist but not project
                if test_poject != False: #if project exist we just have to add the sample
                    sample = str(row[20]) #We create the sample
                    # Load samples_id file and create a samples_id list
                    if args.samples:
                        samples = check_samples(sample,samples,sample_list)
                    else:
                        if sample not in samples: #We check if the samples is already in the sample list
                            samples.append(sample)
                    project["SAMPLES"] = samples #We add the sample to the prohect
                else: #if project does not exist, we have to add the project and the sample
                    project["PROJECT_NAME"] = row[21]
                    projects.append(project)
                    run["PROJECTS"] = projects
                    sample = str(row[20]) #We create the sample
                    # Load samples_id file and create a samples_id list
                    if args.samples:
                        samples = check_samples(sample,samples,sample_list)
                    project["SAMPLES"] = samples #We add the sample to the prohect
            else: #if run does not exist, we have to create the run, the project and the sample
                projects=[]
                samples=[]
                run = {}
                project={}
                run["RUN_NAME"] = row[22] # We create RUN_NAME
                project["PROJECT_NAME"] = row[21] #We create the project
                sample = str(row[20]) #We create the sample
                if args.samples:
                    samples = check_samples(sample,samples,sample_list)
                else:
                    if sample not in samples: #We check if the samples is already in the sample list
                        samples.append(sample)
                project["SAMPLES"] = samples #We add the sample to the prohect
                projects.append(project) #We add project to the run
                run["PROJECTS"] = projects
                runs.append(run) #We add run to run list
            json_data["RUNS"] = runs
        db.close() # disconnect from server

    ## Save json file
    if args.save_json:
        with open(args.save_json + '/' + json_data['SERVICE_NUMBER'] + '_data_template.json', 'w') as json_file:
          json.dump(json_data, json_file, indent=2)
          json_file.write('\n')
    ## Create markdown
    outputText = template.render(json_data)
    file = open(args.output + "/INFRES_"+json_data['SERVICE_NUMBER']+".md","wb")
    file.write(outputText.encode("utf-8"))
    file.close()

    logger.info("Infres has being successfully rendered.")

if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        logger.exception(e)
        raise
