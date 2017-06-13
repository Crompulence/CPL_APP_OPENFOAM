cd ./test_forces_case
    python clean.py
cd -  
rm cpl/coupler_header cpl/map_* &> /dev/null
rm *.dat &> /dev/null
rm *.pyc &> /dev/null
exit 0
