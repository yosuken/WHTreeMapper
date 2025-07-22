
require 'rake'

idx, minseqlen, name, fque, fa, fjsn = ARGV ### fin, fout, fjson
flog = "#{fa}.log"
MinSeqLen = minseqlen.to_i
fwl = open(flog, "w")

# {{{ write_ent(fque, fw, fwl, id, seq, numseq, numex)
def write_ent(fque, fw, fwl, id, seq, numseq, numex)
  seq = seq.join
  len = seq.size
  if len >= MinSeqLen
    numseq += 1 
    fw.puts [">"+id+" ", seq] 

    ### [!!!] space at the end of comment line (">"+id+" ") is VERY IMPORTANT because chunkify parse last underscore + digit (e.g., _123) as abundance
    ### see https://github.com/lczech/gappa/wiki/Subcommand:-chunkify
  else
    numex += 1
  end

  if (numseq + numex) > 0 and (numseq + numex) % 1_000_000 == 0 
    fwl.puts "[#{Time.now.strftime("%Y-%m-%d_%H:%M:%S")}] #{fque}: #{numex + numseq} entries parsed."
    fwl.flush
  end

  return [numseq, numex]
end
# }}}

# {{{ parse_fasta(fque, fw)
def parse_fasta(fque, fw, fwl)
  fwl.puts "[#{Time.now.strftime("%Y-%m-%d_%H:%M:%S")}] start parsing #{fque}" 
  fwl.flush
  numseq, numex = 0, 0
  ids = {}

  ### process each line
  id, seq = "", []
  open(fque){ |fr|
    while l = fr.gets
      if l[0] == ">"
        numseq, numex = write_ent(fque, fw, fwl, id, seq, numseq, numex) if id != ""
        id, seq = "", [] ### re-init

        lab = l[1..-1]
        id = lab.split(/\s+/)[0]

        ### [!!!] id duplication check
        raise("#{Errmsg} sequence id #{id} is found twice. Please ensure that all sequence ids in query files are unique.") if ids[id]
        ids[id] = 1
      else
        seq << l.gsub(/\s+/, "")
      end
    end
    numseq, numex = write_ent(fque, fw, fwl, id, seq, numseq, numex) if id != ""
  }

  return [numseq, numex]
end
# }}}

numseq, numex = 0, 0
open(fa, "w"){ |fw|
  numseq, numex = parse_fasta(fque, fw, fwl)
}

h = { idx: idx.to_i, name: name, numseq: numseq, numtooshortseq: numex, fasta: fa, fjson: fjsn, original: File.absolute_path(fque) }
open(fjsn, "w"){ |fw| fw.puts h.inspect }

fwl.close
