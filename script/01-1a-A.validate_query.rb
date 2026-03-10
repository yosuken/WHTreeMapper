
STDOUT.sync = true; STDERR.sync = true

require 'rake'
require 'json'
require 'digest'

minseqlen, name, fque, fa, fjsn, flst = ARGV
MinSeqLen = minseqlen.to_i
flog = "#{fa}.log"
fwlog = open(flog, "w")

# {{{ def write_ent(fque, fw, fwlog, id, seq, numseq, numex)
def write_ent(fque, fw, fwlog, id, seq, numseq, numex)
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
    fwlog.puts "[#{Time.now.strftime("%Y-%m-%d_%H:%M:%S")}] #{fque}: #{numex + numseq} entries parsed."
    fwlog.flush
  end

  return [numseq, numex]
end
# }}}

# {{{ def parse_fasta(fque, fw, fwlog, fwl)
def parse_fasta(fque, fw, fwlog, fwl)
  fwlog.puts "[#{Time.now.strftime("%Y-%m-%d_%H:%M:%S")}] start parsing #{fque}" 
  fwlog.flush
  numseq, numex = 0, 0
  ids = {}

  ### process each line
  id, seq = "", []

  ### check if fque is gzipped
  is_fque_gz = false
  if fque =~ /\.gz$/ or fque =~ /\.gzip/
    is_fque_gz = true
    require 'zlib'
  end

  fr = is_fque_gz ? Zlib::GzipReader.open(fque) : open(fque)

  while l = fr.gets
    if l[0] == ">"
      numseq, numex = write_ent(fque, fw, fwlog, id, seq, numseq, numex) if id != ""
      id, seq = "", [] ### re-init

      lab = l[1..-1]
      id = lab.split(/\s+/)[0]

      ### [!!!] id duplication check
      raise("#{Errmsg} sequence id #{id} is found twice. Please ensure that all sequence ids in query files are unique.") if ids[id]
      ids[id] = 1
      fwl.puts id
    else
      seq << l.gsub(/\s+/, "")
    end
  end
  numseq, numex = write_ent(fque, fw, fwlog, id, seq, numseq, numex) if id != ""

  fr.close

  return [numseq, numex]
end
# }}}

numseq, numex = 0, 0
fw  = open(fa, "w")
fwl = open(flst, "w")
numseq, numex = parse_fasta(fque, fw, fwlog, fwl)
[fw, fwl, fwlog].each{ |_fw| _fw.close }

### write md5
fasta_MD5 = Digest::MD5.file(fa).to_s
fasta_ori_MD5 = Digest::MD5.file(fque).to_s

h = {
  name: name, numseq: numseq, numtooshortseq: numex, fasta: fa, fjson: fjsn, fasta_ori: File.absolute_path(fque), fasta_MD5: fasta_MD5, fasta_ori_MD5: fasta_ori_MD5
}
open(fjsn, "w"){ |fw| fw.puts h.to_json }
