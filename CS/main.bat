:: 文字コードをUTF-8に変更
chcp 65001

SET /p answer="main.cppを実行しますか(Y/N)?"
if /i {%ANSWER%}=={y} (goto :yes)
if /i {%ANSWER%}=={yes} (goto :yes)
EXIT

:yes
del result\prof\*.prof
del result\vtu\*.vtu
.\build\main.exe 2>result/error.log


cd /Users/takinosenaoto/Downloads/my_research/CS用/input && 
g++ -I/opt/homebrew/include/eigen3 -fopenmp make_input.cpp -o make_input && ./make_input 2>error.log

cd /Users/takinosenaoto/Downloads/my_research/CS用/src && 
g++ -I/opt/homebrew/include/eigen3 -fopenmp main.cpp -o main && ./main 2>error.log
:: .\main.exe >null 2>error.log
./main 2>error.log