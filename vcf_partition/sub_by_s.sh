
l=$1
sub() {
    local result
    result=$(( $1 - $l + 1  ))
    echo "$result"
}

add() {
    local result
    result=$(( $1 + 1 ))
    echo "$result"
}

for p in 1 2 3 4
do
#for l in 1000
#do
	n=$(( $p * $l ))
	sed -n "$(sub $n),$n p" allind > s$p.$l
	bcftools view -S s$p.$l d3.vcf.gz > s$p.$l.vcf

#done
done

#sed -n '1,1000p' physpos > p1_$l
#sed -n '1001,2000p' physpos > p2_$l
#sed -n '2001,3000p' physpos > p3_$l
#sed -n '3001,4000p' physpos > p4_$l
#bcftools filter -r chr10:$(sub $(head -n 1 p1))-$(add $(tail -n 1 p1)) d1.vcf.gz  > p1.vcf
#bcftools filter -r chr10:$(head -n 1 p1_$l)-$(tail -n 1 p1) d1.vcf.gz  > p1.vcf
#bcftools filter -r chr10:$(head -n 1 p2_$l)-$(tail -n 1 p2) d1.vcf.gz  > p2.vcf
#bcftools filter -r chr10:$(head -n 1 p3_$l)-$(tail -n 1 p3) d1.vcf.gz  > p3.vcf
#bcftools filter -r chr10:$(head -n 1 p4)-$(tail -n 1 p4) d1.vcf.gz  > p4.vcf

