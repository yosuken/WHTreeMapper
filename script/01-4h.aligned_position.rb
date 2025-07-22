
odir, fpos, fallaln, frefaln, fasn = ARGV

## parse fpos
lab2pos = {}
IO.readlines(fpos).each.with_index(1){ |l, idx|
  next if l =~ /^#/ ## skip comment lines (begins with "#")

  lab, pos = l.chomp.split("\t")
  lab2pos[lab] = pos.split(",").map{ |i|
    if i.strip !~ /^[0-9]+$/
      raise("\e[1;31mError:\e[0m invalid position value in line #{idx} of #{fpos} - '#{l}'. postion should be positive integer (and comma separeted if multiple values)")
    end
    i.strip.to_i
  }
}

## parse fasn --- assign result (might not exist --> [!!!] in this case, gid2info is empty) 
gid2info = {} ## (name LWR fract aLWR afract taxopath) --> parse fract, taxpath

gid2set  = {} ## consider redundant sequence (ERR315859.153160.1_2_5_1 has redundant mapping result 'ERR315859.153160.1_2_5_1 ;ERR315859.32048204.1_3_6_1')
if File.exist?(fasn)
  IO.readlines(fasn)[1..-1].each{ |l|
    gids, fract, tax = l.chomp.split("\t").values_at(0, 2, 5)

    # gids = gids.split(";").map{ |gid| ## parse ';' separated ids
    #   gid.strip.split(/\s+/)[0] ## exclude additional info
    # }
    gids = gids.split(" ;").map{ |gid| gid.strip }

    gids.each{ |gid|
      gid2info[gid] ||= [] ## store [fract, tax]
      gid2info[gid] << [fract, tax]

      gid2set[gid] = gids
    }
  }

  gid2info.each_key{ |gid|
    a = gid2info[gid]
    best = a.sort_by{ |i| - i[0].to_f }[0] ## extract if fract (fract of LWR) is higher
    gid2info[gid] = best ## only highest likelihood hit
  }
else 
  $stderr.puts "assign result: #{fasn} does not exist. generate result without taxon/clade assignment info."
end

## parse last gid in frefaln
ref_last = ""
IO.read(frefaln).split(/^>/)[1..-1].each{ |ent|
  lab, *seq = ent.split("\n")
  gid = lab.split(/\s+/)[0]
  ref_last = gid ## lastly, the last gid will be stored
}

## parse fallaln
open("#{odir}/aligned_position.tsv", "w"){ |fw|
  fw.puts ["query", lab2pos.keys, "fract", "taxpath"]*"\t"

  flag = 0 ## ref or que

  IO.read(fallaln).split(/^>/)[1..-1].each{ |ent|
    lab, *seq = ent.split("\n")
    gid = lab.split(/\s+/)[0]
    seq = seq.join.gsub(/\s+/, "")

    if flag == 1
      out = [gid]

      ## parse amino acids in specified position
      lab2pos.each{ |lab, pos| out << pos.map{ |i| seq[i-1] }*"," }

      ## add taxon/clade info
      out << (gid2info[gid] ? gid2info[gid] : ["NA", "NA"])

      gid, *info = out

      raise("gid: #{gid} - not found in gid2set[gid].") unless gid2set[gid]
      gid2set[gid].each{ |_gid| ### make multiple lines if member of a set is multiple.
        fw.puts [_gid, info]*"\t"
      }
    end

    flag = 1 if ref_last == gid ## make flag of ref --> que
  }
}
