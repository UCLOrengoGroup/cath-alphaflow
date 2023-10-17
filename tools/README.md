Install Chainsaw

```
cd tools
git clone git@github.com:JudeWells/chainsaw.git
cd chainsaw
python3.8 -m venv venv
. venv/bin/activate
pip install --upgrade pip wheel
pip install -r requirements.txt
```

Install Merizo

```
cd tools
git clone https://github.com/psipred/Merizo.git
cd Merizo
python3.8 -m venv venv
. venv/bin/activate
pip install --upgrade pip wheel
pip install -r requirements.txt
```

Install UniDoc

```
cd tools
wget http://yanglab.nankai.edu.cn/UniDoc/download/UniDoc.tgz
mv UniDoc.tgz UniDoc.tar
tar xvf UniDoc.tar
```
