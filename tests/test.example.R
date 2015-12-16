test.examples <- function()
{
  #checkEquals(6, yguR(2,6))
  #checkEqualsNumeric(6, yguR(2,6))
  #checkIdentical(6, yguR(2,6))
  checkTrue(is.list(seqGen(100,5)))
  checkException(log('a'), 'Unable to take the log() of a string')
}

test.deactivation <- function()
{
  DEACTIVATED('Deactivating this test function')
}