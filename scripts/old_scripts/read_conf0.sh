################### Read the config file ###################
read_conf() {
    local conf=$1;
    tmpFile=$2
    echo "###tmp.txt" > $tmpFile
    while read curline; do
        if [[ $curline != \#* && ! -z $curline ]]; then
            var=`echo $curline | awk -F= '{print $1}'`
            #var=${var// /}
            var=`echo $var | awk '{$1=$1;print}'` ##remove white space at the begin/end

            val=`echo $curline | awk -F= '{print $2}'`
            val=`echo $val | awk -F# '{print $1}'`
            val=`echo $val | awk '{$1=$1;print}'` ##remove white space at the begin/end

            #export $var="$val"
            echo $var="$val" >> $tmpFile
        fi
    done < "$conf"
    source $tmpFile
    
}

