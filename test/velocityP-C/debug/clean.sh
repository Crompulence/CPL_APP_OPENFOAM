cd ./test_vels_case
    python2 clean.py
cd -  
rm cpl/coupler_header cpl/map_* &> /dev/null
rm *.dat &> /dev/null
rm *.pyc &> /dev/null
rm PlyParser* &> /dev/null
exit 0
