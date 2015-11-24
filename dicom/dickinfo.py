__author__ = 'p0054421'
import numpy as np
import os
import os.path
import re
import shutil
import time
from os.path import expanduser

import pydicom

from anonymisator import anonymize

patate = '/home/p0054421/Downloads/testDicom/'


def contained_dirs(dir):
    return filter(os.path.isdir,
                  [os.path.join(dir, f) for f in os.listdir(dir)])


def anon(fname):
    patient_list = contained_dirs(fname)

    for i in patient_list:
        if not os.listdir(i):
            continue
        elif os.path.split(os.path.split(i)[1])[1][:4] == 'AAA_':
            acquis_list = contained_dirs(i)
            for j in acquis_list:
                f = []
                for (dirpath, dirnames, filenames) in os.walk(j):
                    f.extend(filenames)
                    break
                ll = 0
                for k in filenames:
                    image = j + '/' + k
                    ds = pydicom.read_file(image)
                    ID = ds.PatientID
                    new_PN = ID
                    new_ID = 'FLOW_' + os.path.split(i)[1] + '_' + ID
                    anonymize(image, image, new_person_name=new_PN,
                              new_patient_id=new_ID, remove_curves=True, remove_private_tags=True)

                    if ll % 50 == 0:
                        print 'anonymized', ll, 'over', len(filenames)
                    ll += 1

        else:
            continue
    return


def rangator(folder, anonimisation=True):
    aquis_list = contained_dirs(folder)
    tmp = expanduser("~") + '/tmp/'
    if not os.path.exists(tmp):
        os.makedirs(tmp)

    print
    print 'images are not deleting fo\' real but are moved to ', tmp
    print ''
    for i in aquis_list:
        if not os.listdir(i):
            print i, "is empty"
            continue
        elif os.path.split(os.path.split(i)[1])[1][:4] == 'AAA_':
            print i, "already dun"
            continue
        else:
            series_list = contained_dirs(i)
            biggus_dickus = 0
            biggus_dickus_fold = 'biggus_dickus_fold'
            for j in series_list:
                if len(os.listdir(j)) > biggus_dickus:
                    biggus_dickus = len(os.listdir(j))
                    biggus_dickus_fold = j
                elif len(os.listdir(j)) <= biggus_dickus:
                    continue

            print 'Largest folder contains', biggus_dickus, "images."
            print 'Largest folder is', biggus_dickus_fold

            for j in series_list:
                if len(os.listdir(j)) == biggus_dickus:
                    continue
                elif len(os.listdir(j)) < biggus_dickus:
                    print j, ' moving to ', tmp, '\n'
                    # shutil.move(j, tmp)
                    shutil.rmtree(j)
                else:
                    print('wut?')
            # image folder of interrst
            acquis_folder = biggus_dickus_fold
            # print acquis_folder
            patient_folder = os.path.split(os.path.split(acquis_folder)[0])[0]
            im = acquis_folder + '/' + \
                 [f for f in os.listdir(acquis_folder) if os.path.isfile(os.path.join(acquis_folder, f))][0]
            imm = pydicom.read_file(im)
            # print im
            name = imm.PatientName
            initiales = ''
            splitname = re.split(r'[\s,.|^/]+', name)
            lenn = len(splitname) - 1
            for j in xrange(len(splitname)):
                if j < lenn:
                    initiales += splitname[j][:1]
                else:
                    initiales += splitname[j][:3]

            date = imm.AcquisitionDate
            patient_folder_new = patient_folder + '/AAA_' + initiales
            subfolder_name = os.path.split(acquis_folder)[0] + '/AAA_' + initiales + '_' + date

            # subfolder/acquisition folder renaming
            if not os.path.exists(subfolder_name):
                os.rename(acquis_folder, subfolder_name)
            else:
                subfolder_name = subfolder_name + '_doublon-' + str(int(np.array([time.clock()])[0] * 10000))
                os.rename(acquis_folder, subfolder_name)

            # parent folder/patient folder renaming or moving if exist
            print subfolder_name
            if os.path.exists(patient_folder_new):
                print 'exist'
                if os.path.exists(patient_folder_new + '/' + os.path.split(acquis_folder)[1]):
                    subfolder_name2 = subfolder_name + '_doublon-' + str(int(np.array([time.clock()])[0] * 10000))
                    os.rename(subfolder_name, subfolder_name2)
                else:
                    subfolder_name2 = subfolder_name
                shutil.move(subfolder_name2, patient_folder_new)
            else:
                print 'not exist'
                os.makedirs(patient_folder_new)
                shutil.move(subfolder_name, patient_folder_new)
            print 'movint to next acquisition'
            print ''

    if anonimisation:
        anon(patate)

    return


rangator(patate, )
