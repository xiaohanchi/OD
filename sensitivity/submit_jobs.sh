#!/bin/bash

n_jobs=42
n_split=5
jobs_per_round=100
short_jobs=(222)
long_jobs=(179 180 182 183 185 186 188 189 191 192)
# run_jobs=(4 9 14 19 24 29 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 52 57 62 67 72 77 82 87 92 97 102 107 112 117 122 127 132 137 142 147 152 157 162 167 169 170 171 172 173 174 175 176 177 178 179 181 184 187 190 193 196 199 202 205 206 207 244 249 254 259 264 269 274 279 284 289 294 299)

for i in $(seq 1 ${n_jobs}); do
# for i in ${run_jobs[@]}; do

	mkdir -p jobtype$i/results
	mkdir -p jobtype$i/raw_results

	cp main.R main_tmp.R
	perl -pi -e "s/sc00/${i}/g" main_tmp.R
	mv main_tmp.R ./jobtype$i/main.R


	cp runjobs.lsf runjobs_tmp.lsf
	perl -pi -e "s/folder00/jobtype$i/g; s/idx00/$i/g" runjobs_tmp.lsf
	if [[ " ${short_jobs[@]} " =~ " ${i} " ]]; then
		perl -pi -e "s/medium/short/g; s/24:00/3:00/g" runjobs_tmp.lsf
	fi

	if [[ " ${long_jobs[@]} " =~ " ${i} " ]]; then
		perl -pi -e "s/medium/long/g; s/24:00/50:00/g" runjobs_tmp.lsf
	fi
	mv runjobs_tmp.lsf ./jobtype$i/runjobs.lsf

done

echo "DONE Creating Files."

for s in $(seq 1 ${n_split}); do

	start_idx=$(( (s - 1) * jobs_per_round + 1 ))
	end_idx=$(( s * jobs_per_round ))
	
	for i in $(seq 1 ${n_jobs}); do
	# for i in ${run_jobs[@]}; do

		bsub -J "outlier_${i}[${start_idx}-${end_idx}]" < ./jobtype$i/runjobs.lsf

	done
	echo "Round $s submitted"
	sleep 2
done
echo "DONE Submitting Jobs."
