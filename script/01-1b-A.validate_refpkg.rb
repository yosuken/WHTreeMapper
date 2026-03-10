
STDOUT.sync = true; STDERR.sync = true

require 'rake'
require 'json'

rpkg, odir, falnO, ftreO, fhmmO, ftaxO, fposO = ARGV
name = File.basename(rpkg)

### perse LENG of fhmm
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

### make a copy of aligned fasta file with only gene ID (remove description for workaround of witch-ng run issue)
faln = "#{odir}/backbone.mfa"
open(faln, "w"){ |fw|
  IO.read(falnO).split(/^>/)[1..-1].each{ |ent|
    lab, *seq = ent.split("\n")
    gid = lab.split(" ")[0]
    fw.puts ">#{gid}\n#{seq*""}"
  }
}

### make a copy of tree, hmm, taxon.tsv, position.tsv
ftre = "#{odir}/backbone.nwk"
fhmm = "#{odir}/backbone.hmm"
ftax = ftaxO ? "#{odir}/taxon.tsv" : nil
fpos = fposO ? "#{odir}/position.tsv" : nil
open(ftre, "w"){ |fw| fw.puts IO.read(ftreO) } if ftre
open(fhmm, "w"){ |fw| fw.puts IO.read(fhmmO) } if fhmm
open(ftax, "w"){ |fw| fw.puts IO.read(ftaxO) } if ftax
open(fpos, "w"){ |fw| fw.puts IO.read(fposO) } if fpos

### for apples-2 mininmum evolution tree
apdir = "#{odir}/for_apples-2" ; mkdir_p apdir unless Dir.exist?(apdir)
ftreME = "#{apdir}/backbone_min_evo.nwk"
flogME = "#{apdir}/backbone_min_evo.log"
ferrME = "#{apdir}/backbone_min_evo.err"
puts "### Generating minimum evolution distanced tree using FastTree for APPLES-2..."
sh "FastTree -nosupport -nome -noml -log #{flogME} -intree #{ftre} < #{faln} > #{ftreME} 2> #{ferrME}"
puts ""

### for pplacer only when IQTREE tree (workaround for potential taxtastic issue)
ppdir = "#{odir}/for_pplacer" ; mkdir_p ppdir unless Dir.exist?(ppdir)
ftreGM = "#{ppdir}/backbone_gamma.nwk"
flogGM = "#{ppdir}/backbone_gamma.log"
ferrGM = "#{ppdir}/backbone_gamma.err"
puts "### Generating gamma tree using FastTree for pplacer (used only when the input tree is IQTREE)..."
sh "FastTree -nosupport -gamma -nome -mllen -log #{flogGM} -intree #{ftre} < #{faln} > #{ftreGM} 2> #{ferrGM}"
sh "pushd #{ppdir} && taxit create -l backbone -P backbone -t #{File.basename(ftreGM)} -s #{File.basename(flogGM)} -f ../#{File.basename(faln)} --stats-type FastTree && mv backbone/* . && rmdir backbone && popd"
puts ""

h = { name: name, refpkg: rpkg, hmmlen: hmmlen, hmmname: hmmname,
  fhmmO: fhmmO, ftaxO: ftaxO, fposO: fposO, falnO: falnO, ftreO: ftreO,
  fhmm: fhmm, ftax: ftax, fpos: fpos, faln: faln, ftre: ftre,
  ftreME: ftreME, flogME: flogME,
  ftreGM: ftreGM, flogGM: flogGM, ppdir: ppdir,
}

fjsn  = "#{odir}/backbone.json"
open(fjsn, "w"){ |fw| fw.puts h.to_json }
