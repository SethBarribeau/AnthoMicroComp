#!usr/bin/bash

#USAGE: bash ProcessSRA.sh <list of SRR samples> <SRA repo> <Batch Name> <czid metadata> <Czid Project name>

#NOTES
#logging in to CZ ID must be completed before the script is commenced using the below command
# > czid login
# or
# > czid login --persistent

echo "Have you logged into CZID?"
read reply
if [[ $reply == y* ]]; then
	echo "Great. Onwards"
else
	echo "This won't work until you're logged in."
	echo "Please login using czid login"
	exit
fi

tot=$(wc -l $1 | awk '{print$1}')
count=1

echo "Beginning to process batch $3. $tot samples ..."

while read s; do
	echo "$count/$tot: Processing $s ... "
	#download SRA from NCBI
	echo "Fetching $s ..."
	prefetch $s >> $3.log
	#check for success
	check=$(grep -c "'$s' was downloaded successfully" $3.log)
	if [[ $check -eq 1 ]]; then
		#unpack to get to the fastq files
	        echo "Unpacking $s..."
		fasterq-dump $2sra/$s.sra >> $3.log
        	#check if it's paired or single
        	fqcheck=$(ls -lR $s*fastq | wc -l)
		#if it's paired ...
        	if [ $fqcheck -eq 2 ]; then
                	echo "Uploading $s (paired reads) to CZID..."
			#change the naming of the files into a format that CZID can read
                       	mv $s*1.fastq "$s"_R1.fastq
                       	mv $s*2.fastq "$s"_R2.fastq
			# upload to CZid project
			# the output of this process is saved into a log file to look for issues should uploading fail
                       	czid short-read-mngs upload-sample  -p '$5' -s "$s" --metadata-csv $4 $s*R1.fastq $s*R2.fastq >> $3.log
			#the most common issue is czid.org logging out the user whilst uploading. Here we check for this happening
			czcheck=$(tail -n 1 $3.log | grep -c "not authenticated with czid")
			#if this has occurred ...
			if [[ $czcheck -eq 1 ]]; then
				echo "CZID has logged out the user. Halting batch script."
				#log the sample that was half - processed in the Failed output file with the apt reason
				printf "%d\t%s\t%s\n" $count $s "CZID log-in timed out" >> $3_Failed.txt
				#log this failed file, and any that were set to follow it, into another batch file
				#then exit the script
				echo $s > $3.UnProcessed.txt
				sed "0,/$s/d" $1 >> $3.UnProcessed.txt
				exit
			else
				echo "$s successfully uploaded to CZID"
                    		echo "Removing sample to conserve SSD..."
				#remove the fastq files from the working directory
				#please not .sra files will persist in the SRA directory and should be removed after this stage of the process is complete
	                      	rm $s*
			fi
        	else
			#uploading single read fastq files
                	echo "Uploading $s (single read) to CZID..."
                	czid short-read-mngs upload-sample  -p '$5' -s "$s" --metadata-csv $4 $s*.fastq >> $3.log
			czcheck=$(tail -n 1 $3.log | grep -c "not authenticated with czid")
			if [[ $czcheck -eq 1 ]]; then
				echo "CZID has logged out the user. Halting batch script."
				printf "%d\t%s\t%s\n" $count $s "CZID log-in timed out" >> $3_failed.txt
				sed "0,/$s/d" $1 > $3.Unprocessed.txt
				exit
			else
				echo "$s successfully uploaded to CZID"
				echo "Removing sample to conserve SSD ..."
				rm $s*
			fi
        	fi

	else
		#prefetch has failed. Log this and carry on
		echo "$s failed to be retreived from NCBI"
		printf "%d\t%s\t%s\n" $count $s "Prefetch failure" >> $3_Failed.txt
	fi
	#this sample has been fully processed where its been successful.
	#log the count and continue to the next sample
	echo "$s processing complete"
	let count=count+1
done < $1

echo "Processing batch $3 is complete."
