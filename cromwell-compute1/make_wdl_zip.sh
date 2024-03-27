export DIR=.
rm $DIR/workflows.zip
cd $DIR/analysis-wdls/definitions/
zip -r $DIR/workflows.zip .
cd $DIR
