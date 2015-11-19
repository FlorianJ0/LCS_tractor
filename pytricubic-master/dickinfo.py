__author__ = 'p0054421'
import os
import os.path
import re
import time
import shutil
from os.path import expanduser

import numpy as np

import pydicom

patate = '/home/p0054421/Downloads/testDicom/'


# a = glob.glob(patate)
def contained_dirs(dir):
    return filter(os.path.isdir,
                  [os.path.join(dir, f) for f in os.listdir(dir)])


def rangator(folder):
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
            for i in xrange(len(splitname)):
                if i < lenn:
                    initiales += splitname[i][:1]
                    print initiales, splitname[i][:1]
                else:
                    initiales += splitname[i][:3]
                    print initiales, splitname[i][:3]

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
            #

            print 'movint to next acquisition'
            print ''
    return
def anonymisator(fname):

    """
    n [71]: imm.PatientName
    Out[71]: 'H863397'
    In [72]: imm.PatientID
    Out[72]: 'FLOW C-SL 863397'
    """
    kepttags = ['SpecificCharacterSet',
                'ImageType',
                'SOPClassUID',
                'SOPInstanceUID',
                'StudyDate',
                'SeriesDate',
                'AcquisitionDate',
                'ContentDate',
                'StudyTime',
                'SeriesTime',
                'AcquisitionTime',
                'ContentTime',
                'Modality',
                'Manufacturer',
                'InstitutionName',
                'StationName',
                'StudyDescription',
                'SeriesDescription',
                'ManufacturerModelName',
                'DerivationDescription',
                'GroupLength',
                'PatientBirthDate',
                'PatientSex',
                'PatientAge',
                'BodyPartExamined',
                'SliceThickness',
                'KVP',
                'DataCollectionDiameter',
                'ProtocolName',
                'ReconstructionDiameter',
                'DistanceSourcetoDetector',
                'DistanceSourcetoPatientGantry',
                'DetectorTilt',
                'TableHeight',
                'RotationDirection',
                'ExposureTime',
                'XRayTubeCurrent',
                'ExposureTime',
                'FilterType',
                'GeneratorPower',
                'FocalSpot',
                'ConvolutionKernel',
                'PatientPosition',
                'GroupLength',
                'StudyInstanceUID',
                'SeriesInstanceUID',
                'StudyID',
                'SeriesNumber',
                'AcquisitionNumber',
                'InstanceNumber',
                'ImagePositionPatient',
                'ImageOrientationPatient',
                'FrameofReferenceUID',
                'PositionReferenceIndicator',
                'SliceLocation',
                'GroupLength',
                'SamplesperPixel',
                'PhotometricInterpretation',
                'Rows',
                'Columns',
                'PixelSpacing',
                'BitsAllocated',
                'BitsStored',
                'HighBit',
                'WindowCenter',
                'WindowWidth',
                'RescaleIntercept',
                'RescaleSlope',
                'GroupLength']


    patient_list = contained_dirs(fname)

    for i in patient_list:
        # print os.listdir(i)
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
                    # print image
                    ds = pydicom.read_file(image)
                    dllist = ds.dir('')
                    ID = ds.PatientID
                    PN = ds.PatientName
                    for t in dllist:
                        if t in kepttags:
                            continue
                        elif t == 'PatientName':
                            ds.PatientName = ID
                        elif t == 'PatientID':
                            ds.PatientID = 'FLOW_' + os.path.split(i)[1] +'_'+ID
                        else:
                            setattr(ds, t, '')
                    if ll % 50 == 0:
                        print 'anonymized', ll, 'over', len(filenames)
                    ll += 1
                    outname = image
                    ds.save_as(outname)

        else:
            continue
    return





rangator(patate)
# out = '/home/p0054421/Downloads/testDicom/anon.dcm'
anonymisator(patate)