################### Read the config file ###################
read_conf() {
    while read curline; do
        if [[ $curline != \#* && ! -z $curline ]]; then
            var=`echo $curline | awk -F= '{print $1}'`
            #var=${var// /}
            var=`echo $var | awk '{$1=$1;print}'` ##remove white space at the begin/end

            val=`echo $curline | awk -F= '{print $2}'`
            val=`echo $val | awk -F# '{print $1}'`
            val=`echo $val | awk '{$1=$1;print}'` ##remove white space at the begin/end

            export "$var"="$val"
        fi
    done < "$1"
}

