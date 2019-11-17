from context import src
from src import test_module

print('passes ../src import')

x = src.test_module.myclass()
x.hello()

