import urllib2
import json
from os.path import expanduser, isfile

def sendFusiRequenst(seqs):
    """ obtain the fusi score as calculated by Fusi et al's webservice 
    >>> sendFusiRequenst([ "GGGAGGCTGCTTTACCCGCTGTGGGGGCGC"])
    """
    keyFname = expanduser("~/.fusiKey.txt")
    if not isfile(keyFname):
        raise Exception("No ~/.fusiKey.txt file found. Request an API key from azimuth@microsoft.com, write it into this file (single line) and retry")

    api_key = open(keyFname, "r").read().strip()
    paramList = [ (seq, "-1", "-1") for seq in seqs]
                        #"Values": [ [ "GGGAGGCTGCTTTACCCGCTGTGGGGGCGC", "-1", "-1" ] ]
    data =  {

            "Inputs": {

                    "input1":
                    {
                        "ColumnNames": ["sequence", "cutsite", "percentpeptide"],
                        "Values": paramList,
                    },        },
                "GlobalParameters": {
    }
        }

    body = str.encode(json.dumps(data))

    url = 'https://ussouthcentral.services.azureml.net/workspaces/ee5485c1d9814b8d8c647a89db12d4df/services/c24d128abfaf4832abf1e7ef45db4b54/execute?api-version=2.0&details=true'
    headers = {'Content-Type':'application/json', 'Authorization':('Bearer '+ api_key)}

    req = urllib2.Request(url, body, headers)

    try:
        response = urllib2.urlopen(req)

        # If you are using Python 3+, replace urllib2 with urllib.request in the above code:
        # req = urllib.request.Request(url, body, headers) 
        # response = urllib.request.urlopen(req)

        return json.loads(response.read())["Results"]["output2"]["value"]["Values"]

    except urllib2.HTTPError, error:
        print("The request failed with status code: " + str(error.code))

        # Print the headers - they include the requert ID and the timestamp, which are useful for debugging the failure
        print(error.info())

        print(json.loads(error.read())) 

