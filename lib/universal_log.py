#!/usr/bin/python

def write_message(message):
    f = open("/share/pipeline_messages", "a")
    f.write(message)


