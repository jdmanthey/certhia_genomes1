for i in {1..100}; do 
	grep 'NPOP1' bs${i}/1/secondary_contact/secondary_contact.bestlhoods > bs_${i}.lhoods
	for rep in {1..20}; do
		grep -v 'NPOP1' bs${i}/${rep}/secondary_contact/secondary_contact.bestlhoods >> bs_${i}.lhoods
	done
done
