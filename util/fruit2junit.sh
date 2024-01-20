#!/bin/bash
# Fruit2junit.sh by Mireille Schneider and Simon Pinches
# Copied, slightly altered by Daan van Vugt

#######
# Script config
###
# Colors
e="\e[1;31m" # Error color
r="\e[m"     # Reset color

# Parsing
log_s="Start of FRUIT summary:"           # Start of fruit log
log_e="-- End of FRUIT summary"           # End of fruit log
ext_s="-- Start of FRUIT extra taste --"  # Start of extra informations
ext_a="-- End of FRUIT extra taste  --"   # End of extra informations
split='\\o\/'                             # Delimiter for extra informations

# Files
output="fruit2junit.xml"                  # Default output file's name
temp_d=$(mktemp fruit2junit.XXXXXX)       # Temporary file for data extraction

#######
# Function : extract data from file
#  extract {input} {start_cut} {end_cut}
###
extract () {
	found=false
	cat $1 | while read -r line ; do
		if [[ $line == *$2* ]]; then    # Start extraction
			found=true
		elif [[ $line == *$3* ]]; then  # Finish extraction
			found=false
		elif [ $found == true ]; then
			echo -e "$line"
		fi
	done
}

#######
# Define Input / Output
###
# Check the number of arguments
if [ $# == 1 ]; then
	input=$1
elif [ $# == 2 ]; then
	input=$1
	output=$2
else
	echo -e "${e}Invalid arguments !\n"
	echo -e "Usage : \n$0 <input file> [output file]$r"
	exit 1
fi

# Check if the file exist
if [ ! -f $input ]; then
	echo -e "${e}File \"${input}\" not found.$r"
	exit 1
fi

#######
# Parse file
###
# Retrieve the order of tests
order=`awk '/ . : successful assert,   F : failed assert/ {getline; getline; print}' $input`

# Check if input file contain extra data and extract it
extra=false
if [[ $(grep -- "$ext_s" $input) ]]; then
	extra=`extract $input "$ext_s" "$ext_a"`
	nb_ex=`echo -e "$extra" | wc -l`
fi

# Extract data from $input to $temp_d
extract $input "$log_s" "$log_e" > $temp_d

# Init result variables
total=`grep -Po "^Total asserts : (\s*)(\d*)" $temp_d | sed 's/[^0-9][^0-9]*//'`
failures=`grep -Po "^Failed        : (\s*)(\d*)" $temp_d | sed 's/[^0-9][^0-9]*//'`

# Failed message
failed_message=`grep -Po "\[_not_set_\]:(.*)" $temp_d`

# Tests loop variables
l=1
content=''

# Loop through tests
for (( i=0; i<${#order}; i++ )); do
	if [ "${order:$i:1}" == "F" ]; then # Failure
	    data=($(echo -e "$failed_message" | awk "NR==$l" | awk -vRS="]" -vFS="[" '{print $ 2}'))
	    content="$content  <testcase name=\"${data[@]:3}\" classname=\"FRUIT\">\n    <failure message=\"Assertion Failed\" type=\"failure\">Reference value = ${data[1]}\nPresent value = ${data[2]}</failure>\n  </testcase>\n"
	    l=$(($l + 1))
	else                              # Successful
		if [ "$extra" != false ] && [ $i -lt $nb_ex ]; then
			temp=$(echo -e "$extra" | awk "NR==$i+1")
			nm=$(echo "$temp" | sed s/${split}.*//g)
			cl=$(echo "$temp" | sed s/.*${split}//g)
			content="$content  <testcase name=\"$nm\" classname=\"$cl\"></testcase>\n"
		else
			content="$content  <testcase name=\"Test n°$i\" classname=\"FRUIT\"></testcase>\n"
		fi
	fi
done

# Remove temporary file
rm $temp_d 

#######
# Write to output file
###
echo "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" > $output
echo "<testsuite name=\"FRUIT tests\" tests=\"$total\" failures=\"$failures\">" >> $output
echo -e "${content}</testsuite>" >> $output

#######
# End of script
###
echo -e "Status  : Exported to ${output}\nTotal   : ${total}\nFailure : ${failures}"
