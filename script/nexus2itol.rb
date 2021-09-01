
require 'pp'
Line_width = 3

fin, fout = ARGV

def parse_tree(s)
  ancestors = []
  t = {}
  tokens = s.gsub(/\s/, "").split(/([;\(\),:])/)
  # p [s, tokens]

  tokens.each.with_index{|token, i|
    case token
    when "(" ## new branchset
      subt = {}
      t["branchset"] = [subt]
      ancestors << t
      t = subt
    when "," ## another branch
      subt = {}
      ancestors[-1]["branchset"] << subt
      t = subt
    when ")" ## optional name next
      t = ancestors.pop
    when ":" ## optional length next
    else
      x = tokens[i-1]
      if %w|) ( ,|.include?(x)
        t["name"] = token
      elsif x == ":"
        if token =~ /\[\&\!([^=]+)=([^\]]+)\]/ ## parse [&!hoge=fuga]
          t["length"] = $`
          t[$1] = $2
        else
          t["length"] = token
        end
      end
    end
  }

  # pp t
  return t
end

def dfs(t, o)
  name  = t["name"]
  color = t["color"]
  branchset = t["branchset"]

  if branchset
    ## parse 'name' for inner node
    if color
      names = []
      branchset.each{ |subt| names << dfs_color(subt, o) }
      o << [names*"|", color]
    end

    ## dfs
    branchset.each{ |subt| dfs(subt, o) }

  else ## leaf
    raise unless name
    o << [name, color] if color
  end
end

def dfs_color(t, o)
  name  = t["name"]
  color = t["color"]
  branchset = t["branchset"]
  if branchset
    subt = branchset[0] ## only 1 child
    dfs_color(subt, o)
  else ## leaf
    raise unless name
    return name
  end
end

def make_output(t)
  out = []
  out << "TREE_COLORS"
  out << "SEPARATOR TAB"
  out << "DATA"

  ### parse color
  o = []
  dfs(t, o)

  # #NODE_ID TYPE   COLOR   LABEL_OR_STYLE SIZE_FACTOR
  # 915|777  branch #00ff00 normal         3
  o.each{ |a|
    name, col = a
    out << [name, "branch", col, "normal", Line_width]*"\t"
  }

  return out
end

def parse_nexus(fin)
  flag = 0
  IO.readlines(fin).each{ |l|
    if l =~ /^BEGIN TREES;/
      flag = 1
    elsif flag == 1 and l =~ /\s+TREE\s+\S+\s*=\s*(.+)$/ ## REE tree1 = <tree>
      return $1
    end
  }
end

# a = []
# a << "(,,(,));"                               ## no nodes are named
# a << "(A,B,(C,D));"                           ## leaf nodes are named
# a << "(A,B,(C,D)E)F;"                         ## all nodes are named
# a << "(:0.1,:0.2,(:0.3,:0.4):0.5);"           ## all but root node have a distance to parent
# a << "(:0.1,:0.2,(:0.3,:0.4):0.5):0.0;"       ## all have a distance to parent
# a << "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"       ## distances and leaf names (popular)
# a << "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;"     ## distances and all names
# a << "((B:0.2,(C:0.3,D:0.4)E:0.5)A:0.1)F;"    ## a tree rooted on a leaf node (rare)
#
# b = []
# b << "(A:0.1[&!color=#aaa],B:0.2[&!color=#bbb],(C:0.3[&!color=#ccc],D:0.4[&!color=#ddd]):0.5[&!color=#eee]);"       ## distances and leaf names (popular)
# b << "(A:0.1[&!color=#aaa],B:0.2[&!color=#bbb],(C:0.3[&!color=#ccc],D:0.4[&!color=#ddd])E:0.5[&!color=#eee])F;"     ## distances and all names
# b << "((B:0.2[&!color=#bbb],(C:0.3[&!color=#ccc],D:0.4[&!color=#ddd])E:0.5[&!color=#eee])A:0.1[&!color=#aaa])F;"    ## a tree rooted on a leaf node (rare)
#
# c = []
# c << "(((Amphioxus_melanopsin:0.6349[&!color=#82bdfe],(zebrafish_melanospin:0.2363[&!color=#81beff],human_Opn4:0.2788[&!color=#84b9fc])1.000:0.4228[&!color=#83bbfd])0.968:0.2447[&!color=#81beff],((Jumping_spider_Rh1:0.30723[&!color=#81bfff],Drosophila_Rh1:0.51021[&!color=#85b7fb])1.000:0.72259[&!color=#9a8ce5],(scallop_SCOP1:0.48196[&!color=#81bfff],(Octopus_rhodopsin:0.0826[&!color=#81bfff],squid_rhodopsin:0.1978[&!color=#81bfff])1.000:0.33096[&!color=#81bfff])0.217:0.12275[&!color=#81bfff])0.739:0.12532[&!color=#81bfff])1.000:0.4825[&!color=#672266],((AmphiOp1:1.07056[&!color=#bd3fbb],scallop_SCOP2:0.98835[&!color=#919fee])0.979:0.36812[&!color=#9f82e0],((cOpn5L1:0.86603[&!color=#85b7fb],(cOpn5m:0.08719[&!color=#81beff],human_Opn5:0.12913[&!color=#82befe])1.000:0.78119[&!color=#82bcfd])0.990:0.35545[&!color=#87b3f9],((Human_peropsin:0.63606[&!color=#83bcfd],((Human_RGR:1.6367[&!color=#b850c6],Squid_retinoshrome:1.9474[&!color=#8da7f3])0.302:0.14827[&!color=#81beff],Amphi_peropsin:0.62189[&!color=#81bfff])0.053:0.15255[&!color=#81bfff])0.700:0.16144[&!color=#81bfff],Jumping_Spoder_peropsin:0.70132[&!color=#82befe])0.987:0.33392[&!color=#82bcfe])0.786:0.11531[&!color=#81bfff])0.451:0.06687[&!color=#81bfff],(((medaka_TMT-opsin1a:0.65391[&!color=#87b3f9],(Human_encephalopsin:0.34831[&!color=#82befe],Zebrafish_Opn3:0.26497[&!color=#81beff])1.000:0.62825[&!color=#939aec])0.925:0.16924[&!color=#81bfff],(Anopheles_GPROP12:0.76802[&!color=#998fe7],((human_red:0.62203[&!color=#85b8fb],(human_blue:0.58934[&!color=#85b6fa],(bovine_rhodosin:0.05546[&!color=#8aadf6],human_rhodopsin:0.02168[&!color=#81bfff])1.000:0.49708[&!color=#86b4f9])0.969:0.23642[&!color=#81bfff])0.884:0.13476[&!color=#81bfff],(salmon_VA_opsin:0.73756[&!color=#82bcfd],lamprey_parapinopsin:0.64095[&!color=#83bcfd])0.669:0.09346[&!color=#81bfff])1.000:0.48042[&!color=#a080df])0.326:0.05138[&!color=#81bfff])0.976:0.2412[&!color=#89aff7],(sea_anemone_opsin:1.14332[&!color=#bc3fba],(box_jellyfish_opsin:0.77201[&!color=#84b8fc],(hydrozaon_jellyfish_opsin:1.01132[&!color=#83bbfd],hydra_opsin:1.06044[&!color=#86b6fa])0.978:0.45694[&!color=#8da7f3])0.980:0.39269[&!color=#722671])0.913:0.27202[&!color=#000000])0.900:0.1686[&!color=#83bcfd]);"
#
# a.each{ |s|
# b.each{ |s|
# c.each{ |s|
#   t   = parse_tree(s)
#   out = make_output(t)
#   puts out
#   break
# }

### parse nexus
s   = parse_nexus(fin)
t   = parse_tree(s)
out = make_output(t)
if fout
  open(fout, "w"){ |fw| fw.puts out }
else
  puts out
end

