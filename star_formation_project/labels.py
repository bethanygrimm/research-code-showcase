import csv
import json
import math

def label_list():
    '''
    Create a JSON file that stores all the sources and their pertinent data, like class, RA/Dec, etc., read from given txt files
    Output JSON file is located in "./json/labels_orig.json"
    '''
    labels = []

    data6 = {}
    data6['header'] = []

    with open ('table6.txt', 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            data6['header'].append(dict(row))

    # print(len(data6['header'])) #indices 6-430 have actual data
    for i in range(6, 431):
        tstr = ((data6['header'][i]['#Table: J/ApJ/890/130/table6.dat (http://cdsarc.cds.unistra.fr)']))
        dstr = tstr.split("|")
        for j in range(len(dstr)):
            dstr[j] = dstr[j].strip()
            # print(dstr)
        c = dstr[1]
        c = "" if (c=="+") else c
        #Calculate inclination angles here.
        #Only table6 has the data needed for calculating inclination angle
        min = float(dstr[12])
        max = float(dstr[11])
        try:
            inc = 90 - (math.asin(min/max) * 180 / math.pi)
        except ZeroDivisionError:
            inc = 90
        labels.append({"Name": dstr[0].upper(), "RADEC": c, "Class": dstr[6], "dBmaj": dstr[11], "dBmin": dstr[12], "dPa": dstr[13], "Inclination": inc})
    #list - we need name, class, dPa, Lbol, Tbol, e_pa
    # print(labels)

    dataE = {}
    dataE['header'] = []

    with open ('tableE.txt', 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            dataE['header'].append(dict(row))

    # print(len(dataE['header'])) #0-379
    for i in range(380):
        added = False
        tstr = ((dataE['header'][i]['Name\tRA_hr\tRA_min\tRA_s\tDEC_deg\tDEC_arcmin\tDEC_arcsec\tFieldA\tDistFA\tFieldV\tDistFV\tClass\tFlux\teFlux\tPFlux\tRMS\tRmaj\tRmin\tPA\tnorm_name\te_bmaj(arcsec)\te_bmin(arcsec)\te_pa(deg)']))
        dstr = tstr.split('\t')
        name = dstr[0].upper()
        for dict_i in labels:
            if (dict_i["Name"] == name):
                dict_i.update({"e_pa (deg)": dstr[22]})
                added = True
        if (not added):
            labels.append({"Name": dstr[0].upper(), "Class": dstr[11], "e_pa (deg)": dstr[22]})

    data8 = {}
    data8['header'] = []

    with open ('table8.txt', 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            data8['header'].append(dict(row))

    # print(len(data8['header'])) #indices 6-481 have actual data
    for i in range(6, 482):
        added = False
        tstr = ((data8['header'][i]['#Table: J/ApJ/890/130/table8.dat (http://cdsarc.cds.unistra.fr)']))
        dstr = tstr.split("|")
        for j in range(len(dstr)):
            dstr[j] = dstr[j].strip()
            # print(dstr)
        name = dstr[0].upper()
        for dict_i in labels:
            if(dict_i["Name"] == name):
                dict_i.update({"Lbol": dstr[2], "Tbol": dstr[3]})
                added = True
        if (not added):
            labels.append({"Name": dstr[0].upper(), "Class": dstr[4], "Lbol": dstr[2], "Tbol": dstr[3]})

    data = {"Data": labels}
    with open('./json/labels_orig.json', 'w') as f:
        json.dump(data,f)

def label_update(name, pab, par, epab, epar, oab, oar, eoab, eoar):
    '''
    Updates sources from labels_orig.json with position and opening angle information.
    Output JSON file is located in "./json/labels.json"

    Inputs:
        name (str): name of source
        pab (float): outflow position angle, blueshifted lobe
        par (float): outflow position angle, redshifted lobe
        epab (float): error of outflow position angle, blueshifted lobe
        epar (float): error of outflow position angle, redshifted lobe
        oab (float): outflow opening angle, blueshifted lobe
        oar (float): outflow opening angle, redshifted lobe
        eoab (float): error of outflow opening angle, blueshifted lobe
        eoar (float): error of outflow opening angle, redshifted lobe
    '''
    data = json.loads(open('./json/labels_orig.json').read())
    l_list = data["Data"]
    has_dict = False
    for dict_i in l_list:
        if dict_i["Name"] == name:
            has_dict = True
            dict_i.update({"Position angle": {"Blue": pab, "Red": par},
                           "Position angle error": {"Blue": epab, "Red": epar},
                           "Opening angle": {"Blue": oab, "Red": oar},
                           "Opening angle error": {"Blue": eoab, "Red": eoar}})
    if not has_dict:
        l_list.append({"Position angle": {"Blue": pab, "Red": par},
                       "Position angle error": {"Blue": epab, "Red": epar},
                       "Opening angle": {"Blue": oab, "Red": oar},
                       "Opening angle error": {"Blue": eoab, "Red": eoar}})

    data = {"Data": l_list}
    with open('./json/labels.json', 'w') as f:
        json.dump(data,f)

def correct_angle(*args):
    '''
    Given inclination angle of the protostellar disk and opening angles of the protostellar outflows, account for inclination correction,
    and add that as another value.
    Output JSON file is located in "./json/labels.json"

    Inputs:
        *args (str): the name of the source for which opening angle is to be inclination-corrected. Updates the json file immediately
            If no arguments passed, all sources are inclination-corrected
    '''
    data = json.loads(open('labels.json').read())
    l_list = data["Data"]

    #just running correct_angle with no arguments corrects all of them
    all = False
    try:
        name = args[0]
    except IndexError:
        name = ""
        all = True

    if(all):
        for i in l_list:
            try:
                inc = i["Inclination"] * math.pi / 180
                eoab = i["Opening error"]["Blue"] * math.pi / 180
                eoar = i["Opening error"]["Red"] * math.pi / 180
                eocb = math.sin(inc)/((1+(math.tan(eoab/2)*math.sin(inc))**2)*math.cos(eoab/2)*math.cos(eoab/2))
                eocr = math.sin(inc)/((1+(math.tan(eoar/2)*math.sin(inc))**2)*math.cos(eoar/2)*math.cos(eoar/2))
                i.update({"Opening corrected error": {"Blue": eocb, "Red": eocr}})
            except KeyError:
                pass
    else:
        has_dict = False
        for dict_i in l_list:
            if dict_i["Name"] == name:
                has_dict = True
                try:
                    inc = dict_i["Inclination"] * math.pi / 180
                    oab = dict_i["Opening angle"]["Blue"] * math.pi / 180
                    oar = dict_i["Opening angle"]["Red"] * math.pi / 180
                    eoab = dict_i["Opening error"]["Blue"] * math.pi / 180
                    eoar = dict_i["Opening error"]["Red"] * math.pi / 180
                    ocb = 2*math.atan(math.tan((oab/2))*math.sin(inc)) * 180 / math.pi
                    ocr = 2*math.atan(math.tan((oar/2))*math.sin(inc)) * 180 / math.pi
                    eocb = math.sin(inc)/((1+(math.tan(eoab/2)*math.sin(inc))**2)*math.cos(eoab/2)*math.cos(eoab/2))
                    eocr = math.sin(inc)/((1+(math.tan(eoar/2)*math.sin(inc))**2)*math.cos(eoar/2)*math.cos(eoar/2))
                    dict_i.update({"Opening corrected": {"Blue": ocb, "Red": ocr}, "Opening corrected error": {"Blue": eocb, "Red": eocr}})
                except KeyError:
                    pass
        if not has_dict:
            pass

    data = {"Data": l_list}
    with open('labels.json', 'w') as f:
        json.dump(data,f)