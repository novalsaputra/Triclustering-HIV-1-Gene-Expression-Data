delta=$1
lambda=$2

python3 Code/runmethod.py $delta $lambda | tee ket1.txt

git config --global user.email "noval.saputra@sci.ui.ac.id"
git config --global user.name "novalsaputra"
git pull origin master

cd Output
folder=$delta-$lambda
rm -rf $folder
mkdir $folder

cd ..
mv gen.txt Output/$folder/gen.txt
mv kondisi.txt Output/$folder/kondisi.txt
mv waktu.txt Output/$folder/waktu.txt
mv msr.txt Output/$folder/msr.txt
mv tqi.txt Output/$folder/tqi.txt
mv ket.txt Output/$folder/ket.txt

git add .
git commit -m 'add $folder' 
git push origin master

