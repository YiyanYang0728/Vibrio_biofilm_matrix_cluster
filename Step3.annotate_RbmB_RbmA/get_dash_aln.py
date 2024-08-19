#!/usr/bin/env python2
####
#NOTE:
#  This script is meant to get people started using the DASH REST API.
#  If implementing DASH as part of an automated process or in publicly-available
#    software it is also important to follow HTTP query best practices in order to
#    account for common network errors.
####
import json
import urllib2

dash_rest_url = "https://sysimm.org/dash/REST1.0"

def query_dash(url):
    print "Querying %s" % (url)
    #Download response from server
    response = urllib2.urlopen(url)
    raw_response = response.read()
    #Split response into lines
    json_object_lines = raw_response.splitlines()
    #Parse JSON objects
    json_objects = [json.loads(line) for line in json_object_lines]
    #Check JSON objects for errors
    for json_object in json_objects:
        if json_object["StatusCode"] > 0:
            code = json_object["StatusCode"]
            message = json_object["StatusMessage"]
            print "Status (%d): %s" % (code, message)
    return json_objects

#Build a query URL to search for 5 DASH representatives for a sequence
pdbid = "1BHE_A"
search_url = "%s/domains?filter=pdbid=%s&limit=10" % (dash_rest_url, pdbid)

#Query DASH for representatives
dash_reps = query_dash(search_url)

g = open("x", 'w')
for dash_rep in dash_reps:
    #Build a query URL to fetch top 5 alignments by score for search results
    alignment_url = "%s/domain_alignments?filter=pdbid1=%s&order=score=desc&limit=10" % (dash_rest_url, dash_rep["PDBID"])
    #Query DASH for alignments
    alignments = query_dash(alignment_url)
    #Output information from alignments
    for alignment in alignments:
        g.write(">DASH|%s\n%s" % (alignment["ID1"], alignment["PRIMS1"])+"\n")
        g.write(">DASH|%s\n%s" % (alignment["ID2"], alignment["PRIMS2"])+"\n")
g.close()
