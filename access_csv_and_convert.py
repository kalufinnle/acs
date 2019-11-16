#!/usr/bin/env python 

#this code implements a function that accepts a path and returns an iterator so that 
#you can access every row of a CSV file, while also converting the n-th item 
#in every row to your desire output. 
import csv

def convertt(input):
#your desired convert function; here, example shows that convertting the third item in every row to an integer
    int(input)

def dataset(path):
    with open(path, 'rU') as data:
        reader = csv.reader(data)
        for row in reader:
            row[2] = convertt(row[2])
            yield row
            
