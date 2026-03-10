
fref, ftmp1, ftmp2 = ARGV

final = "" ### store gene ID of the last reference
open(fref){ |fr|
  while l = fr.gets
    if l[0] == ">"
      gid = l[1..-1].split(/\s+/)[0]
      final = gid
    end
  end
}
open(ftmp2, "w"){ |fw|
  open(ftmp1){ |fr|
    flag = 0
    _gid = ""
    while l = fr.gets
      if l[0] == ">"
        gid = l[1..-1].split(/\s+/)[0]

        flag = 1 if final == _gid ### if previous gene ID is the last reference, start writing
        _gid = gid ### update _gid
      end
      fw.puts l if flag == 1
    end
  }
}

