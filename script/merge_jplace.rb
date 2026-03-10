
### merge jplace files
### usage: ruby <this scrript> <jplace files>

### [!!! CAUTION]
### This script has strict assumption on input files (e.g.,  all jplace files has the same "fields" field). See below.
### 1. tree field should be identical as string.
### 2. fields field shuold be identical.

fins = ARGV

require 'json'

h = {}
fins.each{ |fin|
  ($stderr.puts "input file (#{fin}) is empty. skip it." ; next) if File.zero?(fin)

  str = IO.read(fin)
  ### '-nan,' --> -0.0
  # jsn = JSON.parse(str, max_nesting: false, allow_nan: true) ### allow_nan: true does not work as expected (ruby 3.4.5)

  # {{{ workaround for -nan in jplace file (curerntly disabled)
  # {{{ example
  # {
  #   "p": [
  #     [ 802, -61463.8, 0.3334, -nan, 6.11352e-06 ],
  #     [ 803, -61463.8, 0.3334, -nan, 6.11352e-06 ],
  #     [ 804, -61463.8, 0.333201, 6.27515e-06, 6.11352e-06 ]
  #   ],
  #   "n": [ "UPI0026439ED4_fm272_to422" ]
  # },
  # }}}

  # if str =~ / -?nan,/
  #   str = str.gsub(/ -?nan,/, " -0.0,")
  #   # [!!!] overwrite the original file (necessary for subsequent runs)
  #   open(fin, "w"){ |fw| fw.puts str }
  # end
  # }}}

  jsn = JSON.parse(str, max_nesting: false)
  # $stderr.puts "parsed: #{fin}"

  if h.size == 0 ### first .jplace file
    h = jsn
  else
    ### Is tree identical?
    raise("tree of #{fin} is different from others. Aborting.") if  h["tree"] != jsn["tree"]
    ### Is fields identical?
    raise("fields of #{fin} is different from others. Even the order shuold be identical. Aborting.") if  h["fields"] != jsn["fields"]
    ### add placements
    h["placements"] += jsn["placements"]
  end
}

puts JSON.pretty_generate(h)
# puts JSON.generate(jsn)
