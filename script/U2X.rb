
fin, fout = ARGV

open(fout, "w"){ |fw|
  open(fin){ |fr|
    while l = fr.gets
      if l[0] == ">"
        fw.puts l
      else
        fw.puts l.gsub("U", "X")
      end
    end
  }
}
