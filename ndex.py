#!/usr/bin/env python

#####################################################
#                                                   #
# ndex.py - Interface to the NDEx server            #
#                                                   #
# (c) 2017, Alberto Riva (ariva@ufl.edu)            #
#           Son Le (nuibang@gmail.com)              #
#           University of Florida                   #
#####################################################

import sys
import requests
from requests_toolbelt.multipart.encoder import MultipartEncoder

class NDEX():
    server = "dev.ndexbio.org"
    username = "ariva"
    password = "b1g:n3ts"

    def __init__(self):
        pass

    def sendGet(self, url, json=False):
        req = "http://" + self.server + "/v2" + url
        r = requests.get(req, auth=(self.username, self.password))
        if json:
            return r.json()
        else:
            return r.text
    
    def sendPut(self, url, body, content_type='text/plain', json=False):
        req = "http://" + self.server + "/v2" + url
        print "Sending..."
        r = requests.put(req, data=body,
                         headers={'Content-Type': content_type})
        print r
        if json:
            return r.json()
        else:
            return r.text

    def sendPost(self, url, body, json=False):
        req = "http://" + self.server + "/v2" + url
        r = requests.post(req, data=body)
        if json:
            return r.json()
        else:
            return r.text

    def createNetwork(self, filename):
        m = MultipartEncoder(fields={'CXNetworkStream': ('filename', open(filename, "rb"), 'text/plain') })
        r = self.sendPut("/network", m, content_type=m.content_type, json=False)

### Main

if __name__ == "__main__":
    N = NDEX()
    #result = N.sendGet("/admin/status", json=True)
    result = N.createNetwork(sys.argv[1])
    print result.status_code
    print result.headers
