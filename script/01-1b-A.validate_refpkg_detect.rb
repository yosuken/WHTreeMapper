
STDOUT.sync = true; STDERR.sync = true

require 'json'

rpkg, odir, falnO, ftreO, fhmmO = ARGV
name = File.basename(rpkg)

### parse LENG and NAME of fhmm
hmmlen = "0"
hmmname = ""
open(fhmmO){ |fr|
  while l = fr.gets
    if l =~ /^LENG\s+(\d+)/
      hmmlen = $1
    elsif l =~ /^NAME\s+(\S+)/
      hmmname = $1
    end
  end
}

### make a copy of aligned fasta file with only gene ID (remove description for workaround)
faln = "#{odir}/backbone.mfa"
open(faln, "w"){ |fw|
  IO.read(falnO).split(/^>/)[1..-1].each{ |ent|
    lab, *seq = ent.split("\n")
    gid = lab.split(" ")[0]
    fw.puts ">#{gid}\n#{seq*""}"
  }
}

### make a copy of tree, hmm
ftre = "#{odir}/backbone.nwk"
fhmm = "#{odir}/backbone.hmm"
open(ftre, "w"){ |fw| fw.puts IO.read(ftreO) } if ftreO
open(fhmm, "w"){ |fw| fw.puts IO.read(fhmmO) } if fhmmO

h = { name: name, refpkg: rpkg, hmmlen: hmmlen, hmmname: hmmname,
  fhmmO: fhmmO, falnO: falnO, ftreO: ftreO,
  fhmm: fhmm, faln: faln, ftre: ftre,
}

fjsn  = "#{odir}/backbone.json"
open(fjsn, "w"){ |fw| fw.puts h.to_json }
