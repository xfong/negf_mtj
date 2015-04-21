#!/bin/bash
#PBS -j oe
#PBS -l host=msee286lnx4.ecn.purdue.edu

export GOPATH=/home/finesse/xfong/proj/min.ecn.purdue.edu/n01/negf_mtj/

go run main.go > /export/min.ecn.purdue.edu/n01/negf_mtj_test_out.txt
