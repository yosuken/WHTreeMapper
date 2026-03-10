
# {{{ example
### [witch-ng input]
# >WP_091443104.1
# ---------------------------------------------------------------M--------------------------------P-LGLYDVVQFAILGAGFALLAYFLYSLTSRDEVS-----AR-YRPSSYAALCLAAVATV-AYLLLYLD------WD-------------------S-------GF--R-L-----E---D--------------------------------------------------------GVY--------------------V-P--------NEE--A--------------R-----T---T-----EGTRYI-DW----SITVPLLTVELLAVCSV--T-G-AAA-R---------RLRSS--------TMAAAFLMIVTGYLGA---Q------V---L----------------D------Q------GR---DR------L--------A----------L---------V-------V---W-----GLI--STAFFA-------YLYVALI--------------GAVRRSLPTM---------G-P-EA---AVSLRNATIV--LLSSFGVYPLVYA----V--------P-V--------------F---------A----------DV----------T--------P-----------A------------W-------F----T-----AMQVGYSAADVVAKIGFGVLVHKVAKLRTAED

### [witch-ng output]
# >WP_091443104.1
# ----------------------------------------------------------------M--------------------------------P-LGLYDVVQFAILGAGFALLAYFLYSLTSRDEVS-----AR-YRPSSYAALCLAAVATV-AYLLLYLD------WD-------------------S-------GF--R-L-----E---D--------------------------------------------------------GVY--------------------V-P--------NEE--A--------------R-----T---T-----EGTRYI-DW----SITVPLLTVELLAVCSV--T-G-AAA-R---------RLRSS--------TMAAAFLMIVTGYLGA---Q------V---L----------------D------Q------GR---DR------L--------A----------L---------V-------V---W-----GLI--STAFFA-------YLYVALI--------------GAVRRSLPTM---------G-P-EA---AVSLRNATIV--LLSSFGVYPLVYA----V--------P-V--------------F---------A----------DV----------T--------P-----------A------------W-------F----T-----AMQVGYSAADVVAKIGFGVLVHKVAKLRTAED-
# }}}

faln, fin, fout = ARGV

N_check = 10000 ### compare between alignments using first N_check sequences

aln0 = []
open(faln){ |fr|
  ### refpkg alignment
  j = 0
  while l = fr.gets
    if l[0] == ">"
      j += 1
      break if j > N_check
    else
      ### sequence should be on a single line
      str = l.strip # read second line (should be the first sequence, not partial)
      (0..str.size-1).each{ |i|
        aln0[i] ||= []
        aln0[i] << str[i]
      }
    end
  end
} 

aln0 = aln0.map{ |x| x.join("") }
L = aln0.size
N = aln0[0].size ### min(N_check, number of sequences in the alignment)

aln1 = []
open(fin){ |fr|
  ### witch-ng alignment
  j = 0
  while l = fr.gets
    if l[0] == ">"
      j += 1
      break if j > N
    else
      ### sequence should be on a single line
      str = l.strip # read second line (should be the first sequence, not partial)
      (0..str.size-1).each{ |i|
        aln1[i] ||= []
        aln1[i] << str[i]
      }
    end
  end
}
aln1 = aln1.map{ |x| x.join("") }
L1 = aln1.size
N1 = aln1[0].size ### min(N_check, number of sequences in the alignment)
raise("Error: unexpected alignment. number of sequences in input alignment (#{N}) and witch-ng output (#{N1}) should be the same") if N != N1
# p ["N", N, "N1", N1, "L", L, "L1", L1]

conv = {} ### aln1 position -> aln0 position
k = 0
(0..L1-1).each{ |i|
  (k..L-1).each{ |j|
    if aln1[i] == aln0[j]
      conv[i] = j
      k = j + 1
      break
    end
  }
}

### now conv is built
### find first position in aln1 that matches aln0[0]

# p conv
if conv.size != L
  raise("Error: only #{conv.size} positions could be mapped out of #{L}, when parsing #{faln} and #{fin}")
end

pos1 = conv.keys
open(fout, "w"){ |fw| 
  open(fin){ |fr|
    while l = fr.gets
      if l[0] == ">"
        fw.puts l
      else
        seq = l.strip.split(//)
        fw.puts seq.values_at(*pos1).join("")
      end
    end
  }
}
