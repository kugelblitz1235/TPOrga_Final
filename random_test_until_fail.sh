#!/bin/bash

while
        ./random_test
        ret=$?
        ((ret == 0))
do
        echo "----------------- Passed tests -----------------"
done
echo "Exit code: $ret"
echo "----------------- Failed tests -----------------"