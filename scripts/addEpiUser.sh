#!/usr/bin/bash
if ! [[ -r "/etc/sssd/sssd.conf" ]]; then
	echo "sudo permission is required!"
	exit 0
fi

while [[ $# -gt 0 ]]
do
key="$1"
case $key in
    -u|--user)
    uName="$2"
    shift # past argument
    ;;
    -i|--uid)
    uID="$2"
    shift # past argument
    ;;
    -t|--tsccUname)
    tscc="$2"
    shift
    ;;
    -h|--help)
    help="TRUE"
    ;;
    *)
    ;;
esac
shift # past argument or value
done
if [ -n "$help" ] || [ -z "$uName" ]
then
  echo "

Adding user to epiglass:
usage:
addEpiUser.sh -u <AD name> [optional]
	-u|--user <value>: provide the AD user name of the user. UCSD authentication is used to grante access to epiglass;
> optional:
	-i|--uid <value>: change the uid on epiglass to match the uid from TSCC, and a folder in scratch folder will be created;
	-t <value>: provide the TSCC user name if it is different from AD user name;
	-h|--help: print this message;
> e.g.
	addEpiUser.sh -u z5ouyang -i 501152 -t ouyang

"
  exit
fi
msg="\nuser $uName will be added to epiglass, a home folder will be created.\n"
if [ ! -z "$uID" ]; then
	msg="${msg}The uid will be set to $uID,"
	if [ -z "$tscc"]; then
		msg="$msg and a scratch folder named $uName will also be created.\n"
	else
		msg="$msg and a scratch folder named $tscc will also be created.\n"
	fi
fi
printf "$msg\n"

printf "Grant $uName to access Epiglass:\n\tsudo sed '/simple_allow_users/s/$/,$uName/' /etc/sssd/sssd.conf\n"
#sudo sed "/simple_allow_users/s/$/,$uName/" /etc/sssd/sssd.conf > /etc/sssd/sssd.conf

homeF="/gpfs/data01/glasslab/home/$uName"
printf "Creat Home Folder for $uName on Epiglass:\n\tsudo mkdir $homeF\n\tsudo chown -R $uName:glass-group $homeF\n"
sudo mkdir $homeF
sudo chown -R $uName:glass-group $homeF

printf "make a short cut to $homeF:\n\tsudo ln -s $homeF /home/$uName\n"
sudo ln -s $homeF /home/$uName

if [ ! -z "$uID" ]; then
	printf "Map epiglass uID to TSCC for $uName:\n\tsudo sss_override user-add $uName -u $uID -g 7300 -h $homeF\n\tsudo systemctl restart sssd\n\tsudo chown -R $uName:glass-group $homeF\n"
	sudo sss_override user-add $uName -u $uID -g 7300 -h $homeF
	sudo systemctl restart sssd
	sudo chown -R $uName:glass-group $homeF
  
	if [ -z "$tscc"]; then
		tscc="$uName"
	fi
	scratchF="/gpfs/data01/scratch/glasslab/$tscc"
	printf "Create a scratch folder named $tscc on epiglass:\n\tsudo mkdir $scratchF\n\tsudo chown -R $uName $scratchF\n\tsudo ln -s $scratchF /data/scratch/$uName\n"
	sudo mkdir $scratchF
	sudo chown -R $uName:glass-group $scratchF
	sudo ln -s $scratchF /data/scratch/$uName
fi

echo "Add $uName successfully"





