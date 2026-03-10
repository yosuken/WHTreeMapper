
fout, *fas = ARGV

AA_Properties = { ## reference: KEGG compound, 2020-05-28
  "K" => {"N" => 1, "S" => 0, "C" => 4, "MW" => 146.1876, }, ## C00047
  "R" => {"N" => 3, "S" => 0, "C" => 4, "MW" => 174.201,  }, ## C00062
  "H" => {"N" => 2, "S" => 0, "C" => 4, "MW" => 155.1546, }, ## C00135
  "D" => {"N" => 0, "S" => 0, "C" => 2, "MW" => 133.1027, }, ## C00049
  "E" => {"N" => 0, "S" => 0, "C" => 3, "MW" => 147.1293, }, ## C00025
  "N" => {"N" => 1, "S" => 0, "C" => 2, "MW" => 132.1179, }, ## C00152
  "Q" => {"N" => 1, "S" => 0, "C" => 3, "MW" => 146.1445, }, ## C00064
  "S" => {"N" => 0, "S" => 0, "C" => 1, "MW" => 105.0926, }, ## C00065
  "T" => {"N" => 0, "S" => 0, "C" => 2, "MW" => 119.1192, }, ## C00188
  "Y" => {"N" => 0, "S" => 0, "C" => 7, "MW" => 181.1885, }, ## C00082
  "A" => {"N" => 0, "S" => 0, "C" => 1, "MW" => 89.0932,  }, ## C00041
  "V" => {"N" => 0, "S" => 0, "C" => 3, "MW" => 117.1463, }, ## C00183
  "L" => {"N" => 0, "S" => 0, "C" => 4, "MW" => 131.1729, }, ## C00123
  "I" => {"N" => 0, "S" => 0, "C" => 4, "MW" => 131.1729, }, ## C00407
  "P" => {"N" => 0, "S" => 0, "C" => 3, "MW" => 115.1305, }, ## C00148
  "F" => {"N" => 0, "S" => 0, "C" => 7, "MW" => 165.1891, }, ## C00079
  "M" => {"N" => 0, "S" => 1, "C" => 3, "MW" => 149.2113, }, ## C00073
  "W" => {"N" => 1, "S" => 0, "C" => 9, "MW" => 204.2252, }, ## C00078
  "G" => {"N" => 0, "S" => 0, "C" => 0, "MW" => 75.0666,  }, ## C00037
  "C" => {"N" => 0, "S" => 1, "C" => 1, "MW" => 121.1582, }, ## C00097
}

header = %w|gene len len_of_std_aa avg_MW N-ARSC C-ARSC S-ARSC K R H D E N Q S T Y A V L I P F M W G C others|
fw = open(fout, "w")
fw.puts header*"\t"

fas.each{ |fa|
  next if !(File.exist?(fa)) or File.zero?(fa)

  IO.read(fa).split(/^>/)[1..-1].each{ |ent|
    lab, *seq = ent.split("\n")
    seq = seq.join.gsub(/\s+/, "").gsub("*", "") ### remove stop codon
    gid = lab.split(/\s+/)[0]

    info  = Hash.new(0)

    ## store info
    info["len"]  = seq.size
    info["gene"] = gid
    count = seq.scan(/./).inject(Hash.new(0)){ |h, i| h[i] += 1; h } 
    # p count
    # p AA_Properties

    AA_Properties.each{ |k, prop|
      info[k]          += count[k]
      info["len_of_std_aa"] += count[k]
      %w|N S C|.each{ |i|
        info["#{i}-ARSC"] += count[k] * prop[i]
      }
      %w|MW|.each{ |i|
        info["avg_#{i}"] += count[k] * prop[i]
      }
    }

    ## store others (i.e., other than 20 standard AAs)
    info["others"] = info["len"] - info["len_of_std_aa"]

    ## divide by len_of_std_aa
    %w|N-ARSC S-ARSC C-ARSC avg_MW|.each{ |i|
      info[i] = "%.6g" % (info[i].to_f / info["len_of_std_aa"])
    }

    fw.puts header.map{ |i| info[i] }*"\t"
  }
}

fw.close
