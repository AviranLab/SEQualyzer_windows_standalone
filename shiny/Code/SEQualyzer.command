# change STR to the path where SEQualyzer located
STR=$(dirname $0)
cd $STR
chmod +x ./config_mac.R
./config_mac.R $STR
