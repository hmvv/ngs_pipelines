import sys
import os
import optparse
import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import scipy.stats as st
import mysql.connector
from mysql.connector import Error
from mysql.connector import errorcode


class Sample:
    def __init__(self, elements):
        self.sampleID = elements[0]
        self.runID = elements[1]
        self.sampleName = elements[2]
        self.assay = elements[3]
        self.instrument = elements[4]


def getConnection(db_host,db,user,passwd):
        try:
            return mysql.connector.connect(host=db_host,database=db, user=user, password=passwd)
        except mysql.connector.Error as error :
            print("Failed to connect database:{}".format(error))
            connection.rollback()

def getSampleInfo(connection,sampleID):
    try:
        cursor = connection.cursor(prepared=True)
        query = "select samples.sampleID, samples.runID, samples.sampleName, assays.assayName, instruments.instrumentName from samples join assays on assays.assayID = samples.assayID join instruments on instruments.instrumentID = samples.instrumentID where samples.sampleID= %s"
        cursor.execute(query, (sampleID,))
        sample = None;
        for row in cursor:
            sample = Sample([el.decode('utf-8') if type(el) is bytearray else el for el in row])
        return sample

    except mysql.connector.Error as error :
        print("Failed to get information from database:{}".format(error))
        connection.rollback()


def updateControlTMBTrend(connection):
    try:
        cursor = connection.cursor(prepared=True)
        query = "select TMBScore from sampleTumorMutationBurden where substring(TMBPair,1,4) = 'COLO'"
        cursor.execute(query)
        scores = []
        for row in cursor:
            scores.append(row[0])
        return scores

    except mysql.connector.Error as error :
        print("Failed to get information from database:{}".format(error))
        connection.rollback()


if len(sys.argv) > 1:
    parser = optparse.OptionParser()
    parser.add_option('-i', '--sampleID')
    parser.add_option('-e', '--environment')
    parser.add_option('-d', '--database_host')
    parser.add_option('-u', '--user')
    parser.add_option('-p', '--password')
    options,args = parser.parse_args()
    sampleID = options.sampleID
    environment = options.environment
    db_host = options.database_host
    db = "ngs_"+environment
    user = options.user
    psswd = options.password

    connection = getConnection(db_host,db,user,psswd)
    sample = getSampleInfo(connection,sampleID)

    if "COLO" in str(sample.sampleName):
        scores = updateControlTMBTrend(connection)
        minscore= np.min(scores) - 1.0
        maxscore= np.max(scores) + 1.0
        plt.ylim(minscore,maxscore)
        sns.lineplot(x=range(1,len(scores)+1) ,y=scores, marker="o",color='r',label="TMB Scores")
        plt.title("TMB Assay Control(COLO829) Trendline, Total:"+str(len(scores)))
        plt.ylabel("TMB Score")
        plt.xlabel("Run Number")
        plt.savefig("/home/environments/"+db+"/assayCommonFiles/tmbAssay/TMB_ControlCOLO829_Scores.png")
        plt.close()

    connection.close()

else:
    print("python tmb_qc_plots.py -h for help")
