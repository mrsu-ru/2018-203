# НЕ НУЖНО УДАЛЯТЬ ФАЙЛЫ ivanovii.* и lab.*

# Введение в выч. методы

// ссылка на удаленный репозиторий 
 git remote add mrsu https://github.com/mrsu-ru/2018-203.git

 // новая ветка (создаем и переключаемся)
 git checkout -b mrsu-master
 
 // затягивает изменения из удаленного репозитория
 git pull mrsu master
 
 // переключаемся на ветку master
 git checkout master
 
 // сливаем ветку mrsu-master в ветку master
 git merge --no-ff mrsu-master
 
 // проталкиваем изменения в удаленный репозиторий origin
 git push origin master
 
 // удаляем ветку mrsu-master
 git branch -d mrsu-master
 
 // лог изменений
 git log --graph --decorate --all
