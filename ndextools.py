#!/usr/bin/env python

#####################################################
#                                                   #
# ndextools.py - Interface to the NDEx server       #
#                                                   #
# (c) 2017, Alberto Riva (ariva@ufl.edu)            #
#           Son Le (nuibang@gmail.com)              #
#           University of Florida                   #
#####################################################

import sys

import ndex.client
from requests.exceptions import HTTPError

class NdexClient():
    ndex = None
    host="http://dev.ndexbio.org"
    username="ariva"
    password="b1g:n3ts"

    def __enter__(self):
        self.ndex = ndex.client.Ndex(host=self.host, username=self.username, password=self.password)
        return self.ndex

    def __exit__(self, a1, a2, a3):
        pass

def searchNetwork(term, account=None):
    with NdexClient() as n:
        nets = n.search_networks(search_string=term, account_name=account)
        return nets

def uploadNetwork(cxfile):
    with NdexClient() as n:
        with open(cxfile, "r") as f:
            url = n.save_cx_stream_as_new_network(f)
        sp = url.rfind("/")
        netid = url[sp+1:]
        summary = n.get_network_summary(netid)
        valid = summary['isValid']
        if valid:
            return netid
        else:
            sys.stdout.write("Error: network `{}' is invalid.\n{}\n".format(netid, summary['errorMessage']))
            return ""

def setNetworkSample(netid, cxfile):
    with open(cxfile, "r") as f:
        cxstring = f.read()

    with NdexClient() as n:
        n.set_network_sample(netid, cxstring)
    return True # n.get_sample_network(netid)

def displaySummary(summ):
    print summ
    for field in ['externalId', 'owner', 'name', 'description', 'nodeCount', 'edgeCount']:
        sys.stdout.write("{}: {}\n".format(field, summ[field]))
    props = summ['properties']
    sys.stdout.write("Properties:\n")
    for prop in props:
        sys.stdout.write("  {predicateString}: {value}\n".format(**prop))

def networkSummary(netid):
    with NdexClient() as n:
        try:
            result = n.get_network_summary(netid)
        except HTTPError as e:
            if e.response.status_code == 404:
                sys.stdout.write("Error: network {} does not exist.\n".format(netid))                
            else:
                sys.stdout.write("Error: {}\n".format(e))
            return None
        displaySummary(result)
        return result

def deleteNetwork(netid):
    with NdexClient() as n:
        n.delete_network(netid)
    return netid

if __name__ == "__main__":
    args = sys.argv[1:]
    if args:
        if args[0] == "upload":
            cxfile = args[1]
            netid = uploadNetwork(cxfile)
            if netid:
                sys.stdout.write("{}\tuploaded\t{}\n".format(netid, cxfile))
        elif args[0] == "search":
            term = args[1]
            if len(args) > 2:
                acc = args[2]
            else:
                acc = None
            result = searchNetwork(term, account=acc)
            nets = result['networks']
            for net in nets:
                sys.stdout.write("{}\t{}\n".format(net['externalId'], net['name']))
        elif args[0] == "summary":
            networkSummary(args[1])
        elif args[0] == "delete":
            sys.stdout.write("{}\tdeleted\n".format(deleteNetwork(args[1])))
        elif args[0] == "sample":
            netid = args[1]
            cxfile = args[2]
            print setNetworkSample(netid, cxfile)
