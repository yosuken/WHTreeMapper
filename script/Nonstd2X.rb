
fin, fout = ARGV

open(fout, "w"){ |fw|
  open(fin){ |fr|
    while l = fr.gets
      if l[0] == ">"
        fw.puts l
      else
        ### U, O, B, Z -> X
        fw.puts l.gsub(/[UOBZ]/, "X")
      end
    end
  }
}
