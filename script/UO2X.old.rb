
fin, fout = ARGV
# k

def write_seq(seq, fw, len)
  _seq = seq.join("").gsub(/[UO]/, "X") 
  _seq = _seq[0, len] if len > 0
  fw.puts _seq 

  len = _seq.size if len == 0
  return len
end

open(fout, "w"){ |fw|
  open(fin){ |fr|
    seq = []
    len = 0
    while l = fr.gets
      if l[0] == ">"
        len = write_seq(seq, fw, len) if seq.size != 0

        fw.puts l
        seq = []
      else
        seq << l.strip
      end
    end
    write_seq(seq, fw, len) if seq.size != 0
  }
}
