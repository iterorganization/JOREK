for i in $(ls ../../testcases/ | grep -v tear_circ_199 ); do 
  if [ -d ../../testcases/$i ]; then
    ./copy.sh tear_circ_199.job $i.job;    
  fi
done
