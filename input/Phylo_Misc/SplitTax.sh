#usr/bin/bash

split -l 50000 rankedlineage_181022.tsv rankedLinSplit_

for f in rankedLinSplit*; do                          
	mv "$f" "$f".tsv
done


