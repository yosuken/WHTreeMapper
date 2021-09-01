
fori, fas, alns, ques, faall, bdir = ARGV

require 'rake'

## fori: alignment given by refpkg
fas  = Dir[fas].sort_by{ |fa| fa.split("/")[-1].split(".")[0].gsub(/^chunk_/, "").to_i } ### fas  = "#{Cnkdir}/#{pkg[:name]}/chunk/chunk_*.fasta"
alns = Dir[alns].sort_by{ |fa| fa.split("/")[-2].gsub(/^chunk_/, "").to_i }              ### alns = "#{Cnkdir}/#{pkg[:name]}/alignment/chunk_*/mafft_add.fa"
ques = ques.split(",")                                                                   ### ques = "#{Resdir}/#{pkg[:name]}/each/query/#{que[:name]/trimmed.fa", ordered same as input query file.

hash2seq = {}
fas.each{ |fa| ## trimmed seq, hashed
  IO.read(fa).split(/^>/)[1..-1].each{ |ent|
    hash, *seq = ent.split("\n")
    seq = seq.join.gsub(/\s+/, "")
    hash2seq[hash] = seq
  }
}

seq2aln = {}
alns.each{ |fa| ## alingned seq, hashed
  IO.read(fa).split(/^>/)[1..-1].each{ |ent|
    hash, *aln = ent.split("\n")
    aln = aln.join.gsub(/\s+/, "")
    seq = hash2seq[hash]
    seq2aln[seq]  = aln
  }
}

header = "" ### alignment given by refpkg
IO.read(fori).split(/^>/)[1..-1].each{ |ent|
  lab, *seq = ent.split("\n")
  header += ">" + lab + "\n" + seq.join.gsub(/\s+/, "") + "\n"
}

fwa = open(faall, "w"); fwa.puts header ## all alignment
ques.each{ |fa| ## trimmed seq, labeled
  name = File.basename(fa).gsub(/\.fa$/, "")
  odir = "#{bdir}/#{name}"; mkdir_p odir
  
  open("#{odir}/aligned.fa", "w"){ |fw| ## each alignment
    fw.puts header

    IO.read(fa).split(/^>/)[1..-1].each{ |ent|
      lab, *seq = ent.split("\n")
      seq = seq.join.gsub(/\s+/, "")
      aln = seq2aln[seq]

      fw.puts [">"+lab, aln]
      fwa.puts [">"+lab, aln]
    }
  }
}
fwa.close
