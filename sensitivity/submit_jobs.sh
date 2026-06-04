#!/bin/bash

n_jobs=224
n_split=5
jobs_per_round=100
short_jobs=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200)
long_jobs=(113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 201 202 203 204 205)
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
