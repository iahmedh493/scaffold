install:
	pip install --upgrade pip &&\
		python -m pip install -r requirements.txt

format:
	black *.py

lint:
	pylint --disable=R,C hello.py

test: 
#pytest
	python -m pytest -vv --cov=hello test_hello.py

