import dicom


def anonymize(filename, output_filename, new_person_name="anonymous",
              new_patient_id="id", remove_curves=True, remove_private_tags=True):
    """Replace data element values to partly anonymize a DICOM file.
    Note: completely anonymizing a DICOM file is very complicated; there
    are many things this example code does not address. USE AT YOUR OWN RISK.
    :param remove_private_tags:
    :param remove_curves:
    :param new_patient_id:
    :param output_filename:
    :param new_person_name:
    :param filename:
    """

    # Define call-back functions for the dataset.walk() function
    def PN_callback(ds, data_element):
        """Called from the dataset "walk" recursive function for all data elements.
        :param data_element:
        :param ds:
        """
        if data_element.VR == "PN":
            data_element.value = new_person_name

    def curves_callback(ds, data_element):
        """Called from the dataset "walk" recursive function for all data elements.
        :param data_element:
        :param ds:
         """
        if data_element.tag.group & 0xFF00 == 0x5000:
            del ds[data_element.tag]

    # Load the current dicom file to 'anonymize'
    dataset = dicom.read_file(filename)

    # Remove patient name and any other person names
    dataset.walk(PN_callback)

    # Change ID
    dataset.PatientID = new_patient_id

    # Remove data elements (should only do so if DICOM type 3 optional)
    # Use general loop so easy to add more later
    # Could also have done: del ds.OtherPatientIDs, etc.
    for name in ['OtherPatientIDs', 'OtherPatientIDsSequence', 'PatientAddress']:
        if name in dataset:
            delattr(dataset, name)

    # # Same as above but for blanking data elements that are type 2.
    # for name in ['PatientBirthDate']:
    #     if name in dataset:
    #         dataset.data_element(name).value = ''

    # Remove private tags if function argument says to do so. Same for curves
    if remove_private_tags:
        dataset.remove_private_tags()
    if remove_curves:
        dataset.walk(curves_callback)

    # write the 'anonymized' DICOM out under the new filename
    # print output_filename
    # print dataset.PatientID
    # print dataset.PatientName
    dataset.save_as(output_filename)
    return
